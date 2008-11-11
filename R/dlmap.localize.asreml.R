`dlmap.localize.asreml` <-
function(input, QTLperChr, ...)
{
  locations <- list()
  results <- list()

  chrPos <- input$chrPos
  mrk.names <- input$mrk.names

  cummrk <- c(0, cumsum(input$n.mrk))

  chrSet <- which(QTLperChr>0)

  # Initialize convergence flag
  results$converge=TRUE

    # Loop over all chromsomes which were significant in the preceding detection step
    for (kk in 1:length(chrSet))
    {
     	no.qtls <- QTLperChr[chrSet[kk]][[1]]
	
  	cat("*******************************************\n")
 	cat("Scanning chromosome ", input$chr.names[chrSet[kk]], " for ", no.qtls, " QTL\n") 
 	cat("*******************************************\n")
	loc1 <- list()
	# 1D scan
	for (mm in 1:no.qtls)
	{
	  # 1D scan for QTL
     	  map.results <- dlmap.map.loc.asreml(input, s.chr=chrSet[kk], chrSet=chrSet, prevLoc=loc1)
	
	  # Update convergence flag
	  if (map.results$converge == FALSE) 	results$converge <- FALSE

	  # Find position of maximum wald score and marker markers
	  # or alternately use the likelihood for localization

	  sel.pos <- which.max(map.results$wald)

    	  mark.l <- as.character(find.flanking(input$genCross, chrSet[kk], chrPos[[paste("chr", chrSet[kk], sep="")]][sel.pos])[,1])
    	  mark.r <- as.character(find.flanking(input$genCross, chrSet[kk], chrPos[[paste("chr", chrSet[kk], sep="")]][sel.pos])[,2])

	  mark.ln <- which(mrk.names==mark.l)-cummrk[chrSet[kk]]
	  mark.rn <- which(mrk.names==mark.r)-cummrk[chrSet[kk]]

	  if (mark.ln==mark.rn) 
	  if ((mark.ln>1)&(mark.rn<input$n.mrk[chrSet[kk]])) {
		mark.ln=mark.ln-1 
		mark.rn=mark.rn+1}

          # For the chromosome being scanned record:
          # 1. pos_chr - actual position of the detected QTL
          # 2. fmrk_chr - flanking markers for the positions of detected QTLs
          locations[[paste("pos_", chrSet[kk], "chr", sep="")]] <- c(locations[[paste("pos_", chrSet[kk], "chr", sep="")]], sel.pos)
          loc1[[paste("pos_", chrSet[kk], "chr", sep="")]] <- c(loc1[[paste("pos_", chrSet[kk], "chr", sep="")]], sel.pos)
          locations[[paste("fmrkL_", chrSet[kk], "chr", sep="")]] <- unique(c(locations[[paste("fmrkL_", chrSet[kk], "chr", sep="")]], mark.ln))
          locations[[paste("fmrkR_", chrSet[kk], "chr", sep="")]] <- unique(c(locations[[paste("fmrkR_", chrSet[kk], "chr", sep="")]], mark.rn))

          loc1[[paste("fmrkL_", chrSet[kk], "chr", sep="")]] <- unique(c(loc1[[paste("fmrkL_", chrSet[kk], "chr", sep="")]], mark.ln))
          loc1[[paste("fmrkR_", chrSet[kk], "chr", sep="")]] <- unique(c(loc1[[paste("fmrkR_", chrSet[kk], "chr", sep="")]], mark.rn))

	  if (mm==1)
	  {
		results$profile[[paste("chr", chrSet[kk],sep="")]] <- rbind(as.vector(chrPos[[paste("chr", chrSet[kk], sep="")]]), map.results$wald)
		rownames(results$profile[[paste("chr", chrSet[kk], sep="")]]) <- c("Position", "Wald")
	  }
  	}

    } # end of loop over significant chromosomes
 
  results$qtls <- locations
  return(results)
}

