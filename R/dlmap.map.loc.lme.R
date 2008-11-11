`dlmap.map.loc.lme` <-
function(input, s.chr, chrSet, prevLoc=NULL)
{
  genCross <- input$genCross
  dfMerged <- input$dfMerged
  fixed <- input$envModel$fixed
  chrPos <- input$chrPos
  mrk.names <- input$mrk.names
  n.pos <- input$n.pos
  n.chr <- input$n.chr
  start.i <- input$start.i

  results <- list()
  f.mrk <- vector()
  f.pos <- vector()
  formula <- list()
  chrRE <- vector()

  cummrk <- c(0, cumsum(input$n.mrk))
  cumpos <- c(0, cumsum(n.pos))
  wald <- vector(length=n.pos[s.chr])

  # Construct vector of chrPos which have already been mapped on othechrSetomosomes
  for (kk in 1:n.chr)
  if (length(prevLoc[[paste("pos_", kk, "chr", sep="")]])!=0)
	f.pos <- c(f.pos, start.i+cumpos[kk]+prevLoc[[paste("pos_", kk, "chr", sep="")]])

  # Construct vector of markers which have already been mapped on othechrSetomosomes
  for (kk in 1:n.chr)
  if (length(prevLoc[[paste("fmrkL_", kk, "chr", sep="")]])>0)
  for (jj in 1:length(prevLoc[[paste("fmrkL_", kk, "chr", sep="")]]))
  	f.mrk <- c(f.mrk, cummrk[kk]+prevLoc[[paste("fmrkL_", kk, "chr", sep="")]][jj]:prevLoc[[paste("fmrkR_", kk, "chr", sep="")]][jj])

  # Set up random effects for markers on each chromosome	
  for (kk in 1:n.chr)
  chrRE[kk] <- paste("pdIdent(~", paste(names(dfMerged)[start.i+sum(n.pos)+setdiff(c((cummrk[kk]+1):cummrk[kk+1]), f.mrk)], collapse="+"), "-1)", sep="")

  # Loop over chrPos on the selected chromosome
  for (jj in 1:n.pos[s.chr])
  {
	# set location and determine flanking markers
	scan.loc <- chrPos[[paste("chr", s.chr, sep="")]][jj]
	mark.l <- as.character(find.flanking(genCross, s.chr, scan.loc)$left)
	mark.r <- as.character(find.flanking(genCross, s.chr, scan.loc)$right)

	mark.ln <- which(mrk.names==mark.l)-cummrk[s.chr]
	mark.rn <- which(mrk.names==mark.r)-cummrk[s.chr]
	wald[jj] <- NA

	# Check that position is not in the same interval as any previously mapped QTL
	if ((!is.element(cummrk[s.chr]+mark.ln, f.mrk))|(!is.element(cummrk[s.chr]+mark.rn, f.mrk)))
	{
  	# Alter chromosome random effects to remove any overlap with fixed effects
 	chrRE[s.chr] <- paste("pdIdent(~", paste(names(dfMerged)[start.i+sum(n.pos)+setdiff((cummrk[s.chr]+1):cummrk[s.chr+1], unique(c(mark.ln, mark.rn, f.mrk)))], collapse="+"), "-1)", sep="")

	# Within each chromosome group, markers are modeled as independent and common variance
	# Only include random effects fochrSetomosomes not being scanned
	if (length(chrSet)>2)
	formula$random <- paste("pdBlocked(list(", paste(chrRE[setdiff(chrSet,s.chr)], collapse=","), "))", sep="")

	if (length(chrSet)==2)
	formula$random <- chrRE[setdiff(chrSet, s.chr)]

	formula$fixed <- paste(as.character(fixed)[2], "~", as.character(fixed)[3], sep="")

	# If there are markers already mapped on othechrSetomosomes, add in as fixed effects 
	if (length(f.pos) >0)
	formula$fixed <- paste(formula$fixed, "+",paste(names(dfMerged)[f.pos], collapse="+"), sep="")
	
	formula$fixed <- paste(formula$fixed, "+", names(dfMerged)[jj+cumpos[s.chr]+start.i], sep="")

	formula$fixgrp <- paste(formula$fixed, "| grp1", sep="")
	formula$fixgrp <- as.formula(formula$fixgrp)
	formula$fixed <- as.formula(formula$fixed)

	gd <- groupedData(formula$fixgrp, data=dfMerged)

	# Fit model - different forms depending on relevant terms
	if (length(chrSet)>1)
	model <- lme(fixed=formula$fixed, random=eval(parse(text=formula$random)), data=gd, control=lmeControl(maxIter=input$maxit), na.action=na.omit)

	if (length(chrSet)==1)
	model <- lme(fixed=formula$fixed, random=~1|grp1, data=gd, control=lmeControl(maxIter=input$maxit), na.action=na.omit)

	# Get output - Wald statistic for position fixed effect
	wald[jj] <- summary(model)$tTable[which(rownames(summary(model)$tTable)==names(dfMerged)[jj+cumpos[s.chr]+start.i]),4]
	} # end of check for distinct intervals
  }

  results$wald <- wald*wald

  return(results)
}

