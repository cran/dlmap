`map.det.asreml` <-
function(input, s.chr, chrSet, prevLoc=NULL, ...)
{
  dfMerged <- input$dfMerged
  envModel <- input$envModel
  n.chr <- length(input$map)
  nphe <- input$nphe

  type <- attr(input, "type") 
  results <- list()
  results$converge <- TRUE
  formula <- list()
  wald <- rep(0, length(input$map[[s.chr]]))

  formula <- envModel

  # Set up random effects for markers on each chromosome	
  for (kk in 1:n.chr)
  formula$group[[paste("g_", kk, "chr", sep="")]] <- setdiff(grep(paste("C", kk, "M", sep=""),colnames(dfMerged)[(nphe+1):ncol(dfMerged)]), match(prevLoc, colnames(dfMerged)[(nphe+1):ncol(dfMerged)])) + nphe

  # Loop over positions on the selected chromosome
  mrkloop <- unlist(lapply(strsplit(names(dfMerged)[setdiff(grep(paste("C", s.chr, "M", sep=""), colnames(dfMerged)[(nphe+1):ncol(dfMerged)]), match(prevLoc, colnames(dfMerged)[(nphe+1):ncol(dfMerged)]))+nphe], "M"), function(x) return(x[2])))

  if (type=="f2")
	mrkloop <- unique(substr(mrkloop, 1, nchar(mrkloop)-1))
  mrkloop <- as.numeric(mrkloop)

  for (jj in 1:length(mrkloop))
  {
	chrnam <- paste("idv(grp(g_", setdiff(chrSet,s.chr), "chr))", sep="")
	formula$random <- paste("~", paste(chrnam, collapse="+"))

	if (!is.null(envModel$random))
	  formula$random <- paste(formula$random, "+", as.character(envModel$random[2]), sep="")
	formula$random <- as.formula(formula$random)

	formula$fixed <- paste(as.character(envModel$fixed)[2], "~", as.character(envModel$fixed[3]), sep="")

	if (length(prevLoc) >0)
	formula$fixed <- paste(formula$fixed, "+", paste(prevLoc, collapse="+"), sep="")
	if (type=="f2")
	formula$fixed <- paste(formula$fixed, "+", paste(paste("C", s.chr, "M", mrkloop[jj], c("D", "A"), sep=""), collapse="+"), sep="") else
	formula$fixed <- paste(formula$fixed, "+", paste("C", s.chr, "M", mrkloop[jj], sep=""), sep="")

	formula$fixed <- as.formula(formula$fixed)

	formula$data <- dfMerged
   	formula$Cfixed <- TRUE
  	formula <- c(formula, ...)
   	formula <- formula[!duplicated(formula)]
	formula <- formula[!sapply(formula, is.null)]
	if (length(chrSet)>1)  
	model <- do.call("asreml", formula)

	if (length(chrSet)==1)
	{
	formula1 <- formula
	formula1$random <- envModel$random
	formula1 <- formula1[!sapply(formula1, is.null)]
	model <- do.call("asreml", formula1)
	}

	if (model$converge==FALSE) 	results$converge <- FALSE

	names <- paste("C", s.chr, "M", mrkloop[jj], c("", "D", "A"), sep="")
	if (model$coefficients$fixed[names(model$coefficients$fixed) %in% names]!=0) 	
	wald[mrkloop[jj]] <- wald.test.asreml(model, list(list(which(model$coefficients$fixed[names(model$coefficients$fixed) %in% names]!=0), "zero")))$zres$zwald
  }

  results$wald <- wald
  return(results)
}

