`dlmap.detect.asreml` <-
function(input, filestem, ...)
{
  chr.names <- input$chr.names
  n.perm <- input$n.perm

  test <- vector()
  loc.temp <- list()
  locations <- list()
  results <- list()
  results$converge=TRUE
  permfile = paste(filestem, ".perm1", sep="")
  logfile = paste(filestem, ".det.log", sep="")

  write("Note: if n.perm=0, p-values and threshold are analytical; if n.perm>0 these are empirical", logfile)
  write("****************************************************", logfile)

  # Fit initial detection step with all chromosomes
  # Output contains p-values for each chromosome variance component
  cat("*******************************************\n")
  cat("* Detection Stage Iteration 1:\n")
  cat("* Testing all chromosomes\n")
  cat("*******************************************\n")

  test <- dlmap.test.asreml(input, chrSet=c(1:input$n.chr), ...)
  if (test$converge == FALSE) 	results$converge <- FALSE
  found <- (min(test$adj.pval) < input$alpha)
  chrSet <- which(test$adj.pval < input$alpha)

  if (n.perm>0)
  { 
    write(paste(c("Perm", chr.names), collapse=" "), permfile)
    write(t(cbind(c(1:n.perm), test$perm.ts)), permfile, append=TRUE)
  }

  # Output observed statistics to logfile
  write(paste("Iteration 1: No. Permutations=",n.perm, sep=""), logfile, append=TRUE)
  write(paste(c("",chr.names), collapse="\t"), logfile, append=TRUE)
  write(paste(c("Obs:",round(test$obs, 4)), collapse="\t"), logfile, append=TRUE)
  write(paste(c("P-val:",round(test$adj.pval, 4)), collapse="\t"), logfile, append=TRUE)
  write(paste(input$alpha*100,"% Genomewide Threshold: ", round(test$thresh, 4), sep=""), logfile, append=TRUE)
  write("Significant chromosomes to be used for scanning/testing:", logfile, append=TRUE)
  write(paste(c("", chr.names[chrSet]), collapse="\t"), logfile, append=TRUE)

  itnum = 1
  while (found==1) 
  {
    itnum=itnum+1
    permfile = paste(filestem, ".perm", itnum, sep="")
    loc.temp <- locations

    markers <- vector(length=length(chrSet))
    for (jj in 1:length(chrSet))
    {
     	no.qtls <- length(locations[[paste("m_", chrSet[jj], "chr", sep="")]])

  	cat("*******************************************\n")
	cat("Marker Selection Iteration ", itnum, " for ", chr.names[chrSet[jj]],"\n")
  	cat("*******************************************\n")
     	map.results <- dlmap.map.det.asreml(input, s.chr=chrSet[jj], chrSet=chrSet, prevLoc=loc.temp, ...)
	
	if (map.results$converge == FALSE) 	results$converge <- FALSE

	sel.mrk <- which.max(map.results$wald)
	markers[jj] <- sel.mrk
	locations[[paste("m_", chrSet[jj], "chr", sep="")]] <- c(locations[[paste("m_", chrSet[jj], "chr", sep="")]], sel.mrk)
        locations[[paste("nq_", chrSet[jj], "chr", sep="")]] <- no.qtls+1

  	if (itnum==2)
	{
	  results$profile[[paste("chr", chrSet[jj],sep="")]] <- rbind(as.vector(input$chrPos[[paste("chr", chrSet[jj], sep="")]]), map.results$wald)
	  rownames(results$profile[[paste("chr", chrSet[jj], sep="")]]) <- c("Position", "Wald")
	}
    } # end of loop over significant chromosomes

    # Output selected markers to logfile
    write(paste(c("Mrk:", markers), collapse="\t"), logfile, append=TRUE)

    # Do not allow more than one QTL per interval on a chromosome
    ex <- vector()
    chrSet.old <- chrSet

    for (ii in 1:length(chrSet))
 	if ((no.qtls+2)>=input$n.mrk[chrSet[ii]])
	  ex <- c(ex, ii)

    chrSet <- setdiff(chrSet, chrSet[ex])
  
    if (length(ex)>0)
    {
	write("*******************************************************", logfile, append=TRUE)
	write("Chromosomes full - One QTL located per interval already", logfile, append=TRUE)
	write(paste("Chromosome(s): ", paste(c("", chr.names[chrSet.old[ex]]), collapse="\t"), " removed from testing set", sep=""), logfile, append=TRUE)
	write("*******************************************************", logfile, append=TRUE)
    }

    if (length(chrSet)>0)
    {
	# Test previous subset of chromosomes for significance 
  	cat("*******************************************\n")
	cat("Testing chromosomes for Iteration ", itnum, "\n")
  	cat("*******************************************\n")
  	test <- dlmap.test.asreml(input, chrSet, prevLoc=locations, ...)

	if (n.perm>0)
	{
  	 write(paste(c("Perm", chr.names[chrSet]), collapse=" "), permfile)
  	 write(t(cbind(c(1:n.perm), test$perm.ts)), permfile, append=TRUE)
	}

	if (test$converge==FALSE) 	results$converge <- FALSE
  	found <- (min(test$adj.pval) < input$alpha)
  	chrSet <- chrSet[which(test$adj.pval < input$alpha)]
    }

  write("***************************************************************", logfile, append=TRUE)
  write(paste("Iteration ", itnum,": No. Permutations=",n.perm, sep=""), logfile, append=TRUE)

  write("Chromosomes from previous iteration: ", logfile, append=TRUE)
  write(paste(c("", chr.names[chrSet.old]), collapse="\t"), logfile, append=TRUE)
  write(paste(c("Obs:",round(test$obs,4)), collapse="\t"), logfile, append=TRUE)
  write(paste(c("P-val:",round(test$adj.pval,4)), collapse="\t"), logfile, append=TRUE)
  write(paste(input$alpha*100,"% Genomewide Threshold: ", round(test$thresh,4), sep=""), logfile, append=TRUE)

  if (length(chrSet)>0)
  {
   write("Significant chromosomes for next round of testing/scanning:", logfile, append=TRUE)
   write(paste(c("",chr.names[chrSet]), collapse="\t"), logfile, append=TRUE)
  }
 }  # End of iterative search

  results$qtls <- locations
  return(results)
}

