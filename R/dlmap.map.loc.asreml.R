`dlmap.map.loc.asreml` <-
function(input, s.chr, chrSet, prevLoc=NULL, ...)
{
  genCross <- input$genCross
  dfMerged <- input$dfMerged
  envModel <- input$envModel
  mrk.names <- input$mrk.names
  chrPos <- input$chrPos
  n.pos <- input$n.pos
  n.chr <- input$n.chr
  start.i <- input$start.i

  results <- list()
  f.mrk <- vector()
  f.pos <- vector()
  formula <- list()
  wald <- vector(length=n.pos[s.chr])

  cummrk <- c(0, cumsum(input$n.mrk))
  cumpos <- c(0, cumsum(n.pos))
 
  for (jj in 1:n.chr)
  if (length(prevLoc[[paste("pos_", jj, "chr", sep="")]])!=0)
	f.pos <- c(f.pos, start.i+cumpos[jj]+prevLoc[[paste("pos_", jj, "chr", sep="")]])

  # Construct vector of markers which have already been mapped on othechrSetomosomes
  for (kk in 1:n.chr)
  if (length(prevLoc[[paste("fmrkL_", kk, "chr", sep="")]])>0)
  for (jj in 1:length(prevLoc[[paste("fmrkL_", kk, "chr", sep="")]]))
  	f.mrk <- c(f.mrk, cummrk[kk]+prevLoc[[paste("fmrkL_", kk, "chr", sep="")]][jj]:prevLoc[[paste("fmrkR_", kk, "chr", sep="")]][jj])

  formula <- envModel

  # Set up random effects for markers on each chromosome	
  for (kk in 1:n.chr)
  formula$group[[paste("g_", kk, "chr", sep="")]] <- start.i+sum(n.pos)+setdiff(c((cummrk[kk]+1):cummrk[kk+1]), f.mrk)

  # Initialize convergence flag for output
  results$converge <- TRUE

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
	chrnam <- paste("idv(grp(g_", setdiff(chrSet,s.chr), "chr))", sep="")
	formula$random <- paste("~", paste(chrnam, collapse="+"))

	# Include spatial/environmental random effects
	if (!is.null(envModel$random))
	  formula$random <- paste(formula$random, "+", as.character(envModel$random[2]), sep="")

	formula$random <- as.formula(formula$random)

	formula$fixed <- paste(as.character(envModel$fixed)[2], "~", as.character(envModel$fixed[3]), sep="")

	if (length(f.pos) >0)
	formula$fixed <- paste(formula$fixed, "+",paste(names(dfMerged)[f.pos], collapse="+"), sep="")
	formula$fixed <- paste(formula$fixed, "+", names(dfMerged)[start.i+jj+cumpos[s.chr]], sep="")
	formula$fixed <- as.formula(formula$fixed)

	formula$data <- input$dfMerged
	formula <- c(formula, ...)
	formula <- formula[!sapply(formula, is.null)]

	if (length(chrSet)>1) model <- do.call("asreml", formula)

	if (length(chrSet)==1)
	{
	formula1 <- formula
	formula1$random <- envModel$random
	formula1 <- formula1[!sapply(formula1, is.null)]
	model <- do.call("asreml", formula1)
	}

	if (model$converge==FALSE) 	results$converge <- FALSE

	wald[jj] <- wald.asreml(model)[which(rownames(wald.asreml(model))==names(dfMerged)[jj+cumpos[s.chr]+start.i]),3]
	} # end of check for distinct intervals
  }

  results$wald <- wald

  return(results)
}

