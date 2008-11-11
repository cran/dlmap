`dlmap.read` <-
function(genfile, phefile, mapfile, phename, filestem)
{
  formula <- list()
  output <- list()

  envData <- read.table(phefile, header=TRUE, check.names=FALSE)
  gen.in <- read.table(genfile, header=TRUE, check.names=FALSE)
  map.in <- read.table(mapfile, header=TRUE, check.names=FALSE)

  if ((ncol(map.in)<2)|(ncol(map.in)>3))
	stop("Mapfile has the wrong number of columns; please examine the documentation for the correct input file format more closely")

  if (any(is.na(match(map.in[,1], colnames(gen.in)[2:ncol(gen.in)]))))
	stop("There are different markers in the map file than in the genotype file") 

  if (any(map.in[,1]!=colnames(gen.in)[2:ncol(gen.in)]))
  {
  	message("Warning: order of markers in map file is different to that in the genotype file. Markers will be rearranged in map order")
	gen.in[,2:ncol(gen.in)] <- gen.in[,match(map.in[,1], colnames(gen.in)[2:ncol(gen.in)])]
	colnames(gen.in)[2:ncol(gen.in)] <- map.in[,1]
  }

  if (is.na(match(phename, colnames(envData))))
	stop("The response ", phename, " is not included in the data; please check names of variables")
 
  idname <- colnames(gen.in)[1]
  if (idname!=colnames(envData)[1]) 
	stop("The first column of the genotype and phenotype files, which should contain the unique identifier, has different names in the two files")

  # check if more than two genotype values 
  alleles <- names(table(unlist(gen.in[,2:ncol(gen.in)])))
  if (length(alleles)>2) stop("More than two genotypes have been input; this is not backcross data")

  message("Genotype ", alleles[1], " will be coded as 1 in the data")
  message("Genotype ", alleles[2], " will be coded as 2 in the data")

  gen.in2 <- apply(gen.in, 2, function(x) return(abbreviate(as.character(x))))
  gen.in[,2:ncol(gen.in)][gen.in2[,2:ncol(gen.in2)]==alleles[1]] <- 0
  gen.in[,2:ncol(gen.in)][gen.in2[,2:ncol(gen.in2)]==alleles[2]] <- 1

  cfile <- paste(filestem, ".chrid.dat", sep="")
  mpfile <- paste(filestem, ".markerpos.txt", sep="")
  mnfile <- paste(filestem, ".mnames.txt", sep="")
  pnfile <- paste(filestem, ".pnames.txt", sep="")
  pfile <- paste(filestem, ".pheno.dat", sep="")
  gfile <- paste(filestem, ".geno.dat", sep="")
  
  # output selected portions of files to create Rqtl object
  write.table(map.in[,2], cfile, row.names=FALSE, col.names=FALSE)
  if (ncol(map.in)==3)
  write.table(map.in[,c(1,3)], mpfile, row.names=FALSE, col.names=FALSE)
  write.table(map.in[,1], mnfile, row.names=FALSE, col.names=FALSE)
  write.table(gen.in[,2:ncol(gen.in)], gfile, row.names=FALSE, col.names=FALSE, na="9")

  write.table(envData[1:nrow(gen.in), which(is.element(colnames(envData), phename))], pfile, na="-", row.names=FALSE, col.names=FALSE)

  write.table(colnames(envData)[which(is.element(colnames(envData), phename))], pnfile, row.names=FALSE, col.names=FALSE)

  if (ncol(map.in)==3)
  rqtlCross <- read.cross(format="gary", genfile=gfile, mapfile=mpfile, phefile=pfile, chridfile=cfile, mnamesfile=mnfile, pnamesfile=pnfile, na.strings="NA")

  if (ncol(map.in)==2)
  rqtlCross <- read.cross(format="gary", genfile=gfile, mapfile=NULL, phefile=pfile, chridfile=cfile, mnamesfile=mnfile, pnamesfile=pnfile)

  colnames(rqtlCross$pheno) <- phename
  rqtlCross$pheno[[idname]] <- gen.in[,1]

  file.remove(gfile, pfile, cfile, pnfile, mnfile, mpfile)

  output$rqtlCross <- rqtlCross
  output$envData <- envData
  output$idname <- idname
  if (ncol(map.in)==2)
  output$estmap <- TRUE

  return(output)
}

