`dlmap.map.det.asreml` <-
function(input, s.chr, chrSet, prevLoc=NULL, ...)
{
  dfMerged <- input$dfMerged
  envModel <- input$envModel
  n.pos <- input$n.pos
  n.chr <- input$n.chr
  n.mrk <- input$n.mrk
  start.i <- input$start.i

  results <- list()
  results$converge <- TRUE
  f.mrk <- vector()
  formula <- list()
  cummrk <- c(0, cumsum(n.mrk))
  wald <- vector(length=n.mrk[s.chr])

  for (kk in 1:n.chr)
  	f.mrk <- c(f.mrk, cummrk[kk]+prevLoc[[paste("m_", kk, "chr", sep="")]])

  formula <- envModel

  # Set up random effects for markers on each chromosome	
  for (kk in 1:n.chr)
  formula$group[[paste("g_", kk, "chr", sep="")]] <- start.i+sum(n.pos)+setdiff(c((cummrk[kk]+1):cummrk[kk+1]), f.mrk)

  # Loop over positions on the selected chromosome
  for (jj in setdiff(c(1:n.mrk[s.chr]), f.mrk-cummrk[s.chr]))
  {
	wald[jj] <- NA

	chrnam <- paste("idv(grp(g_", setdiff(chrSet,s.chr), "chr))", sep="")
	formula$random <- paste("~", paste(chrnam, collapse="+"))

	if (!is.null(envModel$random))
	  formula$random <- paste(formula$random, "+", as.character(envModel$random[2]), sep="")
	formula$random <- as.formula(formula$random)

	formula$fixed <- paste(as.character(envModel$fixed)[2], "~", as.character(envModel$fixed[3]), sep="")

	if (length(f.mrk) >0)
	formula$fixed <- paste(formula$fixed, "+",paste(names(dfMerged)[f.mrk+sum(n.pos)+start.i], collapse="+"), sep="")
	formula$fixed <- paste(formula$fixed, "+", names(dfMerged)[jj+cummrk[s.chr]+sum(n.pos)+start.i], sep="")
	formula$fixed <- as.formula(formula$fixed)

	formula$data <- dfMerged
  	formula <- c(formula, ...)
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

	wald[jj] <- wald.asreml(model)[which(rownames(wald.asreml(model))==names(dfMerged)[jj+cummrk[s.chr]+sum(n.pos)+start.i]),3]
  }

  results$wald <- wald
  return(results)
}

