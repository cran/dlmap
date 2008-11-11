`dlmap.lme` <-
function(genfile="dlgenin.dat", phefile="dlphein.dat", mapfile="dlmapin.dat", phename, fixed=NULL, step=0, fixpos=0, seed=1, maxit=60, alpha=.05, filestem="dl")
{
  if (!require(nlme)) 
	stop("nlme package must be installed to use this function. Please install it from CRAN before proceeding")

  fixed.forma <- NULL
  set.seed(seed)

  if (length(phename)!=1)
  	stop("Can only model one phenotypic trait")

  if (is.null(fixed))
	fixed.forma <- as.formula(paste(phename, "~1", sep=""))

  if (!is.null(fixed)&(length(as.character(fixed))!=2))
  {
	print("Fixed formula has incorrect form: reverting to default of ~1")
	fixed.forma <- as.formula(paste(phename, "~1", sep=""))
  }
  
  if (length(grep("\\|", as.character(fixed))>0))
	stop("Error: Fixed formula cannot contain grouping levels")

  if ((length(fixed)==2)&(as.character(fixed)[1]=="~"))
	fixed.forma <- as.formula(paste(phename, as.character(fixed)[1], as.character(fixed)[2], sep=""))

  data.in <- dlmap.read(genfile, phefile, mapfile, phename, filestem)

  genCross <- data.in$rqtlCross
  envData <- data.in$envData
  idname <- data.in$idname

  input <- dlmap.merge(genCross=genCross, step=step, fixpos=fixpos, nperm=0, alpha=alpha, phename=phename, idname=idname)

  input$envModel <- list()
  input$envModel$fixed <- fixed.forma
  input$dfMerged$grp1 <- 1
  start.i <- input$start.i

  #####
  # Recentering
  #####
  n.mrk <- input$n.mrk
  n.pos <- input$n.pos
  dfMdim <- start.i+sum(n.pos)+sum(n.mrk)

  if (length(which(n.mrk==1))>0)
	message("Warning: Linkage groups ", which(n.mrk==1), " contain only one marker")

   input$dfMerged[, (start.i+1):dfMdim] <- t(t(input$dfMerged[,(start.i+1):dfMdim]) - apply(input$dfMerged[,(start.i+1):dfMdim], 2, function(x) return(mean(x,na.rm=TRUE))))

  message("Beginning detection stage...")
  det.out <- dlmap.detect.lme(input, filestem)

  # List of no. detected QTL from output of detection step
  detqtl <- vector(length=input$n.chr)

  for (i in 1:input$n.chr)
  {
	detqtl[i] <- 0
	if (length(det.out$qtls[[paste("nq_", i, "chr", sep="")]]) > 0)
	detqtl[i] <- det.out$qtls[[paste("nq_",i,"chr", sep="")]]
  }
  message("Detection stage complete. ", sum(detqtl), " total QTL detected")
  loc.out <- NULL
  if (((step>0)|(fixpos>0))&(sum(detqtl)>0))
  {
    message("Beginning localization stage...")
    loc.out <- dlmap.localize.lme(input, QTLperChr=detqtl)
  }
  message("Fitting final model with all QTL")
  cumpos <- c(0, cumsum(n.pos))
  cummrk <- c(0, cumsum(n.mrk))
  final <- list()

  f.pos <- vector()
  c.pos <- vector()

  if (!is.null(loc.out))
  for (ii in 1:input$n.chr)
  if (length(loc.out$qtls[[paste("pos_", ii, "chr", sep="")]])!=0)
  {
        f.pos <- c(f.pos, start.i+cumpos[ii]+loc.out$qtls[[paste("pos_", ii, "chr", sep="")]])
	c.pos <- c(c.pos, input$chrPos[[paste("chr", ii, sep="")]][loc.out$qtls[[paste("pos_", ii, "chr",sep="")]]])
  }
  
  if (is.null(loc.out))
  for (ii in 1:input$n.chr)
  if (length(det.out$qtls[[paste("m_", ii, "chr", sep="")]])!=0)
  {
        f.pos <- c(f.pos, start.i+cummrk[ii]+sum(n.pos)+det.out$qtls[[paste("m_", ii, "chr", sep="")]])
	c.pos <- c(c.pos, find.markerpos(input$genCross, input$mrk.names[det.out$qtls[[paste("m_", ii, "chr",sep="")]]+cummrk[ii]])[,2])
  }
  
  final$fixed <- fixed.forma
  if (length(f.pos)>0)
	final$fixed <- paste(as.character(fixed.forma)[2], "~", as.character(fixed.forma)[3], "+", paste(names(input$dfMerged)[f.pos], collapse="+"), sep="")

  final$fixed <- as.formula(final$fixed)
  mod.new <- lm(final$fixed, data=input$dfMerged, na.action=na.omit)

  if (length(f.pos)>0)
  {
    finalest.new <- vector(length=length(f.pos))
    finalsd.new <- vector(length=length(f.pos))
    for (jj in 1:length(f.pos))
    {
	finalest.new[jj] <- mod.new$coefficients[which(names(mod.new$coefficients)==names(input$dfMerged)[f.pos[jj]])]
	finalsd.new[jj] <- summary(mod.new)$coefficients[which(rownames(summary(mod.new)$coefficients)==names(input$dfMerged)[f.pos[jj]]), 2]
    }
  }

   output <- list()

   output$no.qtl <- sum(detqtl)
   output$final.model <- mod.new
 
   output$cross <- input$genCross
   output$trait <- phename

   if (length(f.pos)>0)
   {
     fmrk <- list()
     fmrk$left <- vector()
     fmrk$right <- vector()
   
     size <- finalest.new
     chr <- rep(input$chr.names, detqtl)
     pos <- c.pos

     zvalue <- finalest.new/finalsd.new
     pvalue <- sapply(zvalue, function(x) 2*(1-pnorm(abs(x))))

     # Wald profile for each chromosome containing QTL
     output$profile <- loc.out$profile
     r.chr <- which(detqtl>0)

     if (!is.null(loc.out))
     for (i in 1:length(r.chr))
     { 
       # Flanking markers for each positioned QTL
       fmrk$left <- c(fmrk$left, input$mrk.names[loc.out$qtls[[paste("fmrkL_", r.chr[i], "chr", sep="")]]+cummrk[r.chr[i]]])
       fmrk$right <- c(fmrk$right, input$mrk.names[loc.out$qtls[[paste("fmrkR_", r.chr[i], "chr", sep="")]]+cummrk[r.chr[i]]])
     } 

     if (is.null(loc.out)){
     for (i in 1:length(r.chr))
     {
	fmrk$left <- c(fmrk$left, input$mrk.names[det.out$qtls[[paste("m_", r.chr[i], "chr", sep="")]]+cummrk[r.chr[i]]])
	fmrk$right <- c(fmrk$right, input$mrk.names[det.out$qtls[[paste("m_", r.chr[i], "chr", sep="")]]+1+cummrk[r.chr[i]]])
     }
     output$profile <- det.out$profile}

     zTable <- as.data.frame(cbind(chr, round(pos,2), fmrk$left, fmrk$right, round(size, 2), round(zvalue, 2), round(pvalue,4)))
     names(zTable) <- c("QTL Chr", "Pos", "Left Flank", "Right Flank", "Size", "Z", "p-value")
     output$zTable <- zTable

   } # End of check for detected QTL
   
   if (sum(detqtl)==0)
 	print("No QTL detected in data. See log file for more details of testing.")

   return(output)
}

