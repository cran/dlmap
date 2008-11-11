`dlmap.merge` <-
function(genCross, envData=NULL, step, fixpos, nperm, alpha, pindex=1, phename, idname, estmap=TRUE) 
{
  if (!is.null(envData))
  if (length(setdiff(genCross$pheno[[idname]], envData[[idname]]))>0)
	message("Warning: Individuals exist with genotype data but no phenotypes and will not be considered in analysis")

  results <- list()
  exp.dat <- fill.geno(genCross, method="argmax")
  n.chr <- nchr(exp.dat)
  npop <-  nind(exp.dat)
  n.mrk <- nmar(exp.dat)

  cummrk <- c(0, cumsum(n.mrk))
  if (estmap)
  {
    message("Map has been replaced by estimated map")
    map <- est.map(exp.dat)            
    exp.dat <- replace.map(exp.dat, map)
  }

  chr.names <- paste("Chr", names(exp.dat$geno), sep="")

  positions <- dlmap.calc.pos(exp.dat, step, fixpos)
  n.pos <- positions$n.pos

  mat.m <- matrix(data=NA, nrow=npop, ncol=sum(n.mrk))

  for (ii in 1:n.chr)
  for (jj in 1:n.mrk[ii])
	mat.m[,jj+cummrk[ii]] <- exp.dat$geno[[ii]]$data[,jj]-1

  if (!is.null(envData))
  dfe <- as.data.frame(envData)
  
  df1 <- as.data.frame(positions$mat)
  names(df1) <- paste('p', c(1:sum(n.pos)), sep="") 
  
  df2 <- as.data.frame(mat.m)
  names(df2) <- paste('m', c(1:sum(n.mrk)), sep="")
  dfMrk <- cbind(exp.dat$pheno[[idname]], df2)
  names(dfMrk)[1] <- idname

  df.anal <- cbind(exp.dat$pheno[[idname]], df1, df2)
  names(df.anal)[1] <- idname

  if (!is.null(envData))
  {
    dfe <- cbind(ord=1:nrow(dfe), dfe)
    df.anal <- merge(dfe, df.anal, by=idname, all.x=TRUE, sort=FALSE)
    df.anal[,(ncol(dfe)+1):ncol(df.anal)] <- apply(df.anal[,(ncol(dfe)+1):ncol(df.anal)], 2, as.numeric)
    df.anal <- df.anal[order(df.anal$ord), -2]
  }

   if (is.null(envData))
   df.anal[[phename]] <- exp.dat$pheno[[phename]]

   results$genCross <- exp.dat
   results$dfMerged <- df.anal
   results$chrPos    <- positions$chrPos
   results$mrk.names <- positions$mnames
   results$chr.names <- chr.names
   results$n.perm <- nperm
   results$n.mrk <- n.mrk
   results$n.pos <- n.pos
   results$n.chr <- n.chr
   results$alpha <- alpha
   results$start.i <- which(names(df.anal)=="p1") - 1
   results$idname <- idname

   if (nperm>0)
   results$dfMrk <- dfMrk

  return(results)
}

