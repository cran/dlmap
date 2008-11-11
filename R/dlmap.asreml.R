`dlmap.asreml` <-
function(genfile="dlgenin.dat", phefile="dlphein.dat", mapfile="dlmapin.dat", phename, baseModel, fixed=NULL, random=NULL, rcov=NULL, sparse=NULL, pedigree, step=0, fixpos=0, seed=1, n.perm=0, alpha=.05, filestem="dl", estmap=TRUE, ...)
{
  if (!require(asreml))
	stop("ASReml must be installed to use this function. Please use dlmap.lme if you do not have a license")

  trace <- paste(filestem, ".trace", sep="")
  ftrace <- file(trace, "w")
  sink(trace, type = "output", append = FALSE)
  on.exit(sink(type = "output"))
  on.exit(close(ftrace), add = TRUE)

  set.seed(seed)

  # Check inputs
  if (!file.exists(genfile))	
	stop("Genotype file ", genfile, " does not exist")
  if (!file.exists(phefile))
	stop("Phenotype file ", phefile, " does not exist")
  if (!file.exists(mapfile))
	stop("Map file ", mapfile, " does not exist")
  if (missing(phename))
	stop("Phenotypic trait name is a required argument and is missing")
  if (step<0)
	stop("Invalid input: step size cannot be < 0")
  if (fixpos<0)
	stop("Invalid input: fixpos cannot be < 0")
  if ((step>0) & (fixpos>0))
	stop("Invalid input: only one of step and fixpos should be > 0")

  if (length(phename)!=1)
	stop("Error: Can only analyze one trait")
 
  fixed.forma <- NULL
  if (is.null(fixed))
  {
	message("Warning: No fixed effects input for spatial model; will
		fit the default model: ", phename,"~1")
	fixed.forma <- as.formula(paste(phename, "~1", sep=""))
  }

  if ((!is.null(fixed))&(length(as.character(fixed))!=2))
  {
	message("Fixed formula has incorrect form: reverting to default of ~1")
	fixed.forma <- as.formula(paste(phename, "~1", sep=""))
  }

  if ((length(fixed)==2)&(as.character(fixed)[1]=="~"))
	fixed.forma <- as.formula(paste(phename, as.character(fixed)[1], as.character(fixed)[2], sep=""))
 

  data.in <- dlmap.read(genfile, phefile, mapfile, phename, filestem)

  if (!is.null(data.in$estmap)) estmap <- TRUE

  genCross <- data.in$rqtlCross
  envData <- data.in$envData
  idname <- data.in$idname

  if (!missing(pedigree))
  if (colnames(pedigree)[1]!=idname)
	stop("Error: The identifier for the pedigree does not match the unique identifier in the genotype file")

  Ainv <- NULL
  if (!missing(pedigree))
  {
	ped.Ainv <- asreml.Ainverse(pedigree)$ginv
	
	Ainv <- list()
	Ainv[[idname]] <- ped.Ainv
  } 


    input <- dlmap.merge(genCross=genCross, envData=envData, step=step, fixpos=fixpos, nperm=n.perm, phename=phename,idname=idname, alpha=alpha, estmap=estmap)
 
  if (missing(baseModel))
  input$envModel <- list(fixed=fixed.forma, random=random, sparse=sparse, rcov=rcov, ginverse=Ainv)
  else input$envModel <- baseModel$call

  start.i <- input$start.i
  input$dfMerged[[idname]] <- as.factor(input$dfMerged[[idname]])

   n.mrk <- input$n.mrk
   n.pos <- input$n.pos
   dimdM <- start.i +sum(n.pos)+sum(n.mrk)

   if (length(which(n.mrk==1))>0)
	message("Warning: Linkage groups ", which(n.mrk==1), " contain only one marker")

   input$dfMerged[, (1+start.i):dimdM] <- t(t(input$dfMerged[,(1+start.i):dimdM]) - apply(input$dfMerged[,(1+start.i):dimdM], 2, function(x) return(mean(x,na.rm=TRUE))))

  if (n.perm>0)
   input$dfMrk[, 2:ncol(input$dfMrk)] <- t(t(input$dfMrk[,2:ncol(input$dfMrk)]) - apply(input$dfMrk[,2:ncol(input$dfMrk)], 2, function(x) return(mean(x,na.rm=TRUE))))

  message("Beginning detection stage...")
  det.out <- dlmap.detect.asreml(input, filestem=filestem, ...)
  detqtl <- vector(length=nchr(input$genCross))

  for (i in 1:nchr(input$genCross))
  {
	detqtl[i] <- 0
	if (length(det.out$qtls[[paste("nq_", i, "chr", sep="")]]) > 0)
	detqtl[i] <- det.out$qtls[[paste("nq_",i,"chr", sep="")]]
  }
   
  message("Detection stage complete. ", sum(detqtl), " total QTL detected")
  loc.out <- NULL
  if ((sum(detqtl)>0)&((step>0)|(fixpos>0)))
  {
    message("Beginning localization stage...")
    loc.out <- dlmap.localize.asreml(input, QTLperChr=detqtl, ...)
  }

  message("Fitting final model with all QTL")
  cumpos <- c(0, cumsum(n.pos))
  cummrk <- c(0, cumsum(n.mrk))

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
  
  final.fixed <- fixed.forma
  if (length(f.pos)>0)
	final.fixed <- paste(as.character(fixed.forma)[2], "~", as.character(fixed.forma)[3], "+", paste(names(input$dfMerged)[f.pos], collapse="+"), sep="")

  final.fixed <- as.formula(final.fixed)
  final <- list(fixed=final.fixed, random=input$envModel$random, sparse=input$envModel$sparse, rcov=input$envModel$rcov, ginverse=Ainv, data=input$dfMerged, ...)
  final <- final[!sapply(final, is.null)]
  mod.fin <- do.call("asreml", final)
 
  # get out size estimates from final model output. next statement is incorrect
  finalest.new <- vector(length=length(f.pos))
  finalsd.new <- vector(length=length(f.pos))
  for (jj in 1:length(f.pos))
  {
	finalest.new[jj] <- mod.fin$coefficients$fixed[which(names(mod.fin$coefficients$fixed)==names(input$dfMerged)[f.pos[jj]])]
	finalsd.new[jj] <- (mod.fin$vcoeff$fixed[which(names(mod.fin$coefficients$fixed)==names(input$dfMerged)[f.pos[jj]])]*mod.fin$sigma2)^.5
  }	

   output <- list()
   output$no.qtl <- sum(detqtl)
   output$final.model <- mod.fin
 
   if (length(f.pos)>0)
   {
     fmrk <- list()
     fmrk$left <- vector()
     fmrk$right <- vector()
     size <- finalest.new
     zvalue <- finalest.new/finalsd.new
     pvalue <- sapply(zvalue, function(x) 2*(1-pnorm(abs(x))))
     chr <- rep(input$chr.names, detqtl)
     pos <- c.pos

     # Wald profile for each chromosome containing QTL
     output$profile <- loc.out$profile
     r.chr <- which(detqtl>0)

     if (!is.null(loc.out))
     for (i in 1:length(r.chr))
     { 
       fmrk$left <- c(fmrk$left, input$mrk.names[loc.out$qtls[[paste("fmrkL_", r.chr[i], "chr", sep="")]]+cummrk[r.chr[i]]])
       fmrk$right <- c(fmrk$right, input$mrk.names[loc.out$qtls[[paste("fmrkR_", r.chr[i], "chr", sep="")]]+cummrk[r.chr[i]]])
     } 

     if (is.null(loc.out)) {
     for (i in 1:length(r.chr))
     {
	fmrk$left <- c(fmrk$left, input$mrk.names[det.out$qtls[[paste("m_", r.chr[i], "chr", sep="")]]+cummrk[r.chr[i]]])
	fmrk$right <- c(fmrk$right, input$mrk.names[det.out$qtls[[paste("m_", r.chr[i], "chr", sep="")]]+1+cummrk[r.chr[i]]])
     }
     output$profile <- det.out$profile}

   zTable <- as.data.frame(cbind(chr, round(pos,2), fmrk$left, fmrk$right, round(size,2), round(zvalue,2), round(pvalue,4)))

   names(zTable) <- c("QTL Chr", "Pos", "Left Flank", "Right Flank", "Size", "Z", "p-value")

   output$zTable <- zTable

   } # End of check for detected QTL
   
   if (sum(detqtl)==0)
 	message("No QTL detected in data. See log file for more details of testing.")
   output$cross <- input$genCross
   output$trait <- phename
 
   return(output)
}

