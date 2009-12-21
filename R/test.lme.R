`test.lme` <-
function(input, chrSet, prevLoc=NULL)
{
  dfMerged <- input$dfMerged
  fixed <- input$envModel$fixed
  formula <- list() 
  random.forma <- list() 
  results <- list()
  f.mrk <- vector()
  n.chrSet <- length(chrSet)
  null.forma <- list()
  null.ll <- vector(length=n.chrSet)
  chrRE <- vector()

  LRTStats <- vector() 

  # Construct vector of already mapped markers (f.mrk) 
  formula$fixed <- paste(as.character(fixed)[2], as.character(fixed)[1], as.character(fixed)[3], sep="")

  # Include fixed effects for all markers which have already been mapped
  if (length(prevLoc)>0)
  formula$fixed <- paste(formula$fixed, "+", paste(prevLoc, collapse="+"))

  formula$fixgrp <- paste(formula$fixed, "|grp1", sep="")
  formula$fixed <- as.formula(formula$fixed)
  formula$fixgrp <- as.formula(formula$fixgrp)

  gd <- groupedData(formula$fixgrp, data=dfMerged)

  # Random effects for all markers on a chromosome excluding those which
  # enter the model as fixed effects
  for (ii in 1:length(input$map))
   	chrRE[ii] <- paste("pdIdent(~", paste(setdiff(names(dfMerged)[grep(paste("C", ii, "M", sep=""), names(dfMerged))], prevLoc), collapse="+"), "-1)", sep="")

  formula$random <- paste("pdBlocked(list(", paste(chrRE[chrSet], collapse=","), "))", sep="")

  if (length(chrSet)==1)
  formula$random <- chrRE[chrSet]

  full <- lme(fixed=formula$fixed, random=eval(parse(text=formula$random)), data=gd, control=lmeControl(maxIter=input$maxit), na.action=na.omit)

  full.ll <- full$logLik

  # If there is only one chromosome in the subset, compare a full model to the model
  # with no random effects
  if (n.chrSet==1)
    null.forma[[1]] <- lme(fixed=formula$fixed, data=gd, control=lmeControl(maxIter=input$maxit), na.action=na.omit)
 
  # Otherwise, compare the full model to leave-one-VC-out models, removing 
  # each chromosome effect one at a time
  if (n.chrSet>1)
  for (cc in 1:n.chrSet)
  {
    random.forma[[cc]] <- paste("pdBlocked(list(", paste(chrRE[setdiff(chrSet, chrSet[cc])], collapse=","), "))", sep="")

    if (n.chrSet==2)
    random.forma[[cc]] <- chrRE[setdiff(chrSet, chrSet[cc])]
  
    # Fit the null model, where we omit the specified chromosome random effect
    null.forma[[cc]] <- lme(fixed=formula$fixed, random=eval(parse(text=random.forma[[cc]])), data=gd, control=lmeControl(maxIter=input$maxit), na.action=na.omit)

    null.ll[cc] <- null.forma[[cc]]$logLik
  }
 
  LRTStats <- 2*(full.ll-null.ll)
  results$obs <- LRTStats
  results$raw.pval <- sapply(LRTStats, pvfx)
  results$adj.pval <- results$raw.pval*n.chrSet
  results$thresh <- qchibar(input$alpha)

  return(results)
}

