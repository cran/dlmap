`test.asreml` <-
function(input, chrSet, prevLoc=NULL, ...)
{
  dfMerged <- input$dfMerged
  envModel <- input$envModel
  n.perm <- input$nperm
  n.chr <- length(input$map)
  nphe <- input$nphe
  idname <- input$idname

  formula <- envModel 
  random.forma <- list() 
  results <- list()
  null.forma <- list()
  n.chrSet <- length(chrSet)
  null.ll <- vector(length=n.chrSet)

  npop <- ngen(input)

  # Create permutation matrices  
  perm.test <- matrix(nrow=n.perm+1, ncol=n.chrSet)
  maxp <- vector(length=n.perm+1)
  perm.mat <- matrix(nrow=npop, ncol=n.perm+1)
  perm.mat[,1] <- c(1:npop)

  if (n.perm>0)
  for (kk in 2:(n.perm+1))
	perm.mat[,kk] <- sample(npop)

  results$converge <- TRUE

  # Fit full model, with all chromosomes having separate VCs
    formula$fixed <- paste(as.character(envModel$fixed)[2], "~",as.character(envModel$fixed[3]), sep="")

  # Include fixed effects for all markers which have already been mapped
  if (length(prevLoc)>0)
  formula$fixed <- paste(formula$fixed, "+", paste(prevLoc, collapse="+"), sep="")
  formula$fixed <- as.formula(formula$fixed)

  # Random effects for all markers on a chromosome excluding those which
  # enter the model as fixed effects
  for (ii in 1:n.chr)
    formula$group[[paste("g_", ii, "chr", sep="")]] <- nphe+setdiff(grep(paste("C", ii, "M", sep=""), names(dfMerged)[(nphe+1):ncol(dfMerged)]), match(prevLoc, names(dfMerged)[(nphe+1):ncol(dfMerged)]))

  # Random effects for each chromosome in selected subset
  # Markers are modelled as independent and same variance within chromosomes
  chrnam <- paste("idv(grp(g_", chrSet, "chr))", sep="")
  formula$random <- paste("~", paste(chrnam, collapse="+"))

  if (length(envModel$random)>0)
	formula$random <- paste(formula$random, "+", as.character(envModel$random[2]), sep="")

  formula$random <- as.formula(formula$random)

  formula$dump.model <- TRUE
  formula$data <- dfMerged
  formula <- c(formula, ...)
  formula <- formula[!duplicated(formula)]
  formula <- formula[!sapply(formula, is.null)]
 
  full <- do.call("asreml", formula)

  if (n.chrSet==1)
  {
     form.null <- formula
     form.null$random <- envModel$random
     form.null <- form.null[!sapply(form.null, is.null)]
     null.forma[[1]] <- do.call("asreml", form.null)
  }
 
  if (n.chrSet>1)
  for (cc in 1:n.chrSet)
  {
 	# fit model leaving out each chromosome to test VC
	chrnam <- paste("idv(grp(g_", setdiff(chrSet, chrSet[cc]), "chr))", sep="")
	rndf <- paste("~", paste(chrnam, collapse="+"), sep="")
	if (!is.null(envModel$random))
	rndf <- paste(rndf, "+", as.character(envModel$random[2]), sep="")

	rndf <- as.formula(rndf)
  	form.null <- formula	
	form.null$random <- rndf
	null.forma[[cc]] <- do.call("asreml", form.null)
  }

  # Vector of observed test statistics from LRTs
  if (n.perm==0)
  {
	run <- asreml(model=full) 
	full.ll <- run$loglik
	if (run$converge==FALSE) results$converge <- FALSE

	for (cc in 1:n.chrSet)
	{
	  run <- asreml(model=null.forma[[cc]])
	  null.ll[cc] <- run$loglik
	  if (run$converge==FALSE) results$converge <- FALSE
 	}
  	perm.test[1,] <- 2*(full.ll-null.ll)
	results$obs <- perm.test[1,]

	results$raw.pval <- sapply(perm.test[1,], pvfx)
	results$adj.pval <- sapply(results$raw.pval, function(x) return(min(x*n.chrSet,1)))
	results$thresh <- qchibar(input$alpha/n.chrSet)
  }

  if (n.perm>0)
  {
    dfMrk <- input$dfMrk
    namesrnd <- setdiff(names(dfMrk)[2:ncol(dfMrk)], prevLoc)

    for (ii in 1:(n.perm+1))
    {
        df2 <- cbind(dfMrk[,1],dfMrk[perm.mat[,ii],2:ncol(dfMrk)])
        names(df2)[1] <- names(dfMrk)[1]

        df4 <- dfMerged[, which(!(names(dfMerged)%in%namesrnd))]
        df3 <- merge(df4, df2, by=idname, all.x=TRUE, sort=FALSE)

        # replace data in model for random marker effects
        index <- match(namesrnd, names(full$data))
        index <- index[!is.na(index)]
        index2 <- match(names(full$data)[index], names(df3))
        full$data[,index] <- df3[, index2]

	# run full model
	run <- asreml(model=full)
	full.ll <- run$loglik
        if (run$converge==FALSE)
	  results$converge <- FALSE

	for (cc in 1:n.chrSet)
	{
	   # replace data for random marker effects for each of the null models
          index <- match(namesrnd, names(null.forma[[cc]]$data))
          index <- index[!is.na(index)]
          index2 <- match(names(null.forma[[cc]]$data)[index], names(df3))
          null.forma[[cc]]$data[, index] <- df3[, index2]
          run <- asreml(model=null.forma[[cc]])
          if (run$converge==FALSE) results$converge <- FALSE
                null.ll[cc] <- run$loglik
	}

	perm.test[ii,] <- 2*(full.ll-null.ll)
    }	# end of loops over permutations

  	# For each permutation, store the maximum (over chromosomes) LRT
	maxp <- apply(perm.test, 1, max)

	# Permutation threshold is the (1-alpha) percentile of max values
  	results$thresh <- sort(maxp[2:(n.perm+1)])[floor((1-input$alpha)*n.perm)]
	results$raw.pval <- sapply(perm.test[1,], pvfx)
  	results$adj.pval <- sapply(perm.test[1,], function(x) sum(x<=maxp[2:(n.perm+1)])/n.perm)

	results$obs <- perm.test[1,]
	results$perm.ts <- perm.test

  } # end of check whether n.perm>0

  return(results)
}

