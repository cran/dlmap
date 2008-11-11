`dlmap.calc.pos` <-
function(exp.dat, step, fixpos)
{
  n.chr <- nchr(exp.dat)
  npop <- nind(exp.dat)
  output <- list()
  mat.scan <- list()
  n.pos <- vector(length=n.chr)
  mrk.names <- vector()
  chrPos <- list()
  cumpos <- vector(length=n.chr+1)
  n.mrk <- nmar(exp.dat)
  cummrk <- c(0, cumsum(n.mrk))

  if (fixpos==0)
  {
      geno.p <- calc.genoprob(exp.dat, step=step)
      cumpos[1] <- 0

      for (ii in 1:n.chr)
      {
        n.pos[ii] <- ncol(geno.p$geno[[ii]]$prob)
        cumpos[ii+1] <- n.pos[ii]+cumpos[ii]
        mat.scan[[paste("chr", ii, sep="")]] <- matrix(data=NA, nrow=npop, ncol=n.pos[ii])
	mat.scan[[paste("chr", ii, sep="")]][,] <- geno.p$geno[[ii]]$prob[,,2]

 	mrk.names <- c(mrk.names, colnames(exp.dat$geno[[ii]]$data))
    	chrPos[[paste("chr",ii, sep="")]] <- attr(geno.p$geno[[ii]]$prob, "map")
      }
  }
 
  if ((step==0)&(fixpos>0))
  {
	mrk.names <- vector(length=sum(n.mrk))
	n.pos <- (fixpos+1)*(n.mrk-1)+1
	cumpos <- c(0, cumsum(n.pos))
  	for (ii in 1:n.chr)
	{
	  tempchross <- subset(exp.dat, ii)
	  tempmap <- pull.map(tempchross)
  	  cummap <- tempmap[[1]][1:n.mrk[ii]]
	  intlgth <- tempmap[[1]][2:n.mrk[ii]]-tempmap[[1]][1:(n.mrk[ii]-1)]
	  intstep <- intlgth/(fixpos+1)
	  mat.scan[[paste("chr", ii, sep="")]] <- matrix(data=NA, nrow=npop, ncol=n.pos[ii])
 	  mrk.names[(cummrk[ii]+1):cummrk[ii+1]] <- paste(colnames(exp.dat$geno[[ii]]$data), sep="")
	  chrPos[[paste("chr", ii, sep="")]] <- vector(length=n.pos[ii])

	  for (jj in 1:(n.mrk[ii]-1))
	  {
		tempchross2 <- tempchross
		if (jj<(n.mrk[ii]-1))
		tempchross2 <- drop.markers(tempchross, names(tempchross$geno[[1]]$map[c((jj+2):n.mrk[ii])]))
		tempchross3 <- tempchross2
		if (jj>1)
		tempchross3 <- drop.markers(tempchross2, names(tempchross2$geno[[1]]$map[c(1:(jj-1))]))

		chrPos[[paste("chr", ii, sep="")]][((fixpos+1)*(jj-1)+2):((fixpos+1)*jj)] <- intstep[jj]*c(1:(fixpos))+cummap[jj]
	
		tempgp <- calc.genoprob(tempchross3, step=intstep[jj])
	        mat.scan[[paste("chr", ii, sep="")]][,(1+(fixpos+1)*(jj-1)):((fixpos+1)*jj)] <- tempgp$geno[[1]]$prob[,1:(fixpos+1),2]

	  } 
	mat.scan[[paste("chr", ii, sep="")]][,(fixpos+1)*(n.mrk[ii]-1)+1] <- tempgp$geno[[1]]$prob[,ncol(tempgp$geno[[1]]$prob),2]

	chrPos[[paste("chr", ii, sep="")]][which((c(1:n.pos[ii])-1)%%(fixpos+1)==0)] <- cummap
  	} 
  } 
 
  output$mnames <- mrk.names
  output$mat <- mat.scan
  output$chrPos <- chrPos
  output$n.pos <- n.pos

  return(output)
}

