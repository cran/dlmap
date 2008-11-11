`dlmap.map.det.lme` <-
function(input, s.chr, chrSet, prevLoc=NULL)
{
  dfMerged <- input$dfMerged
  fixed <- input$envModel$fixed
  n.pos <- input$n.pos
  n.mrk <- input$n.mrk
  n.chr <- input$n.chr
  start.i <- input$start.i

  results <- list()
  f.mrk <- vector()
  formula <- list()
  chrRE <- vector()
  wald <- vector(length=n.mrk[s.chr])

  cummrk <- c(0, cumsum(n.mrk))

  for (kk in 1:n.chr)
  	f.mrk <- c(f.mrk, cummrk[kk]+prevLoc[[paste("m_", kk, "chr", sep="")]])

  # Set up random effects for markers on each chromosome	
  for (kk in 1:n.chr)
  chrRE[kk] <- paste("pdIdent(~", paste(names(dfMerged)[start.i+sum(n.pos)+setdiff((cummrk[kk]+1):cummrk[kk+1], f.mrk)], collapse="+"), "-1)", sep="")

  for (jj in setdiff(c(1:n.mrk[s.chr]), f.mrk-cummrk[s.chr]))
  {
	wald[jj] <- NA

  	# Alter chromosome random effects to remove any overlap with fixed effects
  	chrRE[s.chr] <- paste("pdIdent(~", paste(names(dfMerged)[start.i+sum(n.pos)+setdiff((cummrk[s.chr]+1):cummrk[s.chr+1], c(jj+cummrk[s.chr], f.mrk))], collapse="+"), "-1)", sep="")

	# Only include random effects for chromosomes not being scanned
	if (length(chrSet)>2)	
	formula$random <- paste("pdBlocked(list(", paste(chrRE[setdiff(chrSet, s.chr)], collapse=","), "))", sep="")

	if (length(chrSet)==2)
	formula$random <- chrRE[setdiff(chrSet,s.chr)]

	formula$fixed <- paste(as.character(fixed)[2], "~", as.character(fixed)[3], sep="")

	# If there are markers already mapped on other chromosomes, add in as fixed effects 
	if (length(f.mrk) >0)
	formula$fixed <- paste(formula$fixed, "+",paste(names(dfMerged)[f.mrk+sum(n.pos)+start.i], collapse="+"), sep="")

	formula$fixed <- paste(formula$fixed, "+", names(dfMerged)[jj+cummrk[s.chr]+sum(n.pos)+start.i], sep="")

	formula$fixgrp <- paste(formula$fixed, "| grp1", sep="")
	formula$fixed <- as.formula(formula$fixed)
	formula$fixgrp <- as.formula(formula$fixgrp)

	gd <- groupedData(formula$fixgrp, data=dfMerged)

	# Fit model - different forms depending on relevant terms
	if (length(chrSet)>1)
	model <- lme(fixed=formula$fixed, random=eval(parse(text=formula$random)), data=gd, control=lmeControl(maxIter=input$maxit), na.action=na.omit)

	if (length(chrSet)==1)
	model <- lme(fixed=formula$fixed, random=~1|grp1, data=gd, control=lmeControl(maxIter=input$maxit), na.action=na.omit)

	# Get output - Wald statistic for position fixed effect
	wald[jj] <- summary(model)$tTable[which(rownames(summary(model)$tTable)==names(dfMerged)[jj+cummrk[s.chr]+sum(n.pos)+start.i]),4]
  }

  results$wald <- wald*wald

  return(results)
}

