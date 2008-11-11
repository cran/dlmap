`dlmap.convert.cross` <-
function(format=c("csv", "csvr", "csvs", "csvsr", "mm", "qtx", "qtlcart", "gary", "karl", "rqtl"), obj, pheobj, idname="ID", dir="", genoutfile="dlgenin.dat", pheoutfile="dlphein.dat", mapoutfile="dlmapin.dat", file, genfile, mapfile, phefile, chridfile, mnamesfile, pnamesfile, na.strings=c("-", "NA"), genotypes=c("A", "H", "B", "D", "C"), alleles=c("A", "B"), ...)
{
  if (missing(format)) stop("format is a required argument and is missing")

  if ((format=="csv") | (format=="csvr"))
  inCross <- read.cross(format=format, dir=dir, file=file, na.strings=na.strings, genotypes=genotypes, alleles=alleles, ...)

  if ((format=="csvs") | (format=="csvsr"))
  inCross <- read.cross(format=format, dir=dir, genfile=genfile, phefile=phefile, na.strings=na.strings, genotypes=genotypes, alleles=alleles, ...)

  if (format=="mm")
  inCross <- read.cross(format=format, dir=dir, file=file, mapfile=mapfile, ...)

  if (format=="qtx")
  inCross <- read.cross(format=format, dir=dir, file=file, alleles=alleles, ...)

  if (format=="qtlcart")
  inCross <- read.cross(format=format, dir=dir, file=file, mapfile=mapfile, alleles=alleles, ...)

  if (format=="gary")
  inCross <- read.cross(format=format, dir=dir, genfile=genfile, mnamesfile=mnamesfile, chridfile=chridfile, phefile=phefile, pnamesfile=pnamesfile, mapfile=mapfile, na.strings=na.strings, alleles=alleles, ...)

  if (format=="karl")
  inCross <- read.cross(format=format, dir=dir, genfile=genfile, mapfile=mapfile, phefile=phefile, alleles=alleles, ...)

  if (format=="rqtl")
  {
	if (missing(obj)) stop("Must input a cross object if format is rqtl")
	if (!inherits(obj,"cross")) stop("Input object does not have class cross")
	inCross <- obj
  }

  if (is.na(match(idname, names(inCross$pheno))))
	stop("ID variable not found in marker data")

  if (!missing(pheobj)) {
  if (is.na(match(idname, names(pheobj))))
	stop("ID variable not found in supplemental environmental object")
  if (length(intersect(names(table(inCross$pheno[[paste(idname)]])), names(table(pheobj[[paste(idname)]]))))==0)
	stop("ID variables in marker and environmental data do not coincide")
  if (is.na(match(idname, names(pheobj))))
	  stop("Error: ID variable is not in environmental dataset")

  idIndex <- match(idname, names(pheobj))
  }

  n.mrk <- sum(nmar(inCross))
  n.pop <- nind(inCross)
  n.phe <- nphe(inCross)
  n.chr <- nchr(inCross)

  chrid <- NULL
  for (i in 1:n.chr) {
        chrname <- names(inCross$geno[i])
        chrid <- c(chrid, rep(chrname, nmar(inCross)[i]))
  }
  markpos <- NULL
  for (i in 1:n.chr) markpos <- c(markpos, inCross$geno[[i]]$map)
  mnames <- names(markpos)
  geno <- NULL
  for (i in 1:n.chr) geno <- cbind(geno, inCross$geno[[i]]$data)
  geno <- geno - 1

  genout <- as.data.frame(matrix(nrow=n.pop, ncol=(n.mrk+1)))
  pheout <- as.data.frame(matrix(nrow=n.pop, ncol=(n.phe)))
  if (!missing(pheobj))
  pheout <- as.data.frame(matrix(nrow=nrow(pheobj), ncol=ncol(pheobj)))
  mapout <- as.data.frame(matrix(nrow=n.mrk, ncol=3))

  mapout[,1] <- mnames
  mapout[,2] <- chrid
  mapout[,3] <- markpos
  colnames(mapout) <- c("MrkID", "Chr", "Pos")

  genout[,1] <- inCross$pheno[[idname]]
  genout[,2:ncol(genout)] <- geno
  colnames(genout) <- c(idname, mnames)

  if (missing(pheobj)) {
  pheout[,1] <- inCross$pheno[[idname]]
  pheout[,2:ncol(pheout)] <- inCross$pheno[,which(names(inCross$pheno)!=idname)]
  colnames(pheout) <- c(idname, setdiff(names(inCross$pheno), idname))
  }
  else { 
  pheout[,1] <- pheobj[,idIndex]
  pheout[,2:ncol(pheout)] <- pheobj[,setdiff(1:ncol(pheobj), idIndex)]
  colnames(pheout) <- c(idname, names(pheobj)[setdiff(1:ncol(pheobj), idIndex)])  }

  write.table(genout, genoutfile, row.names=FALSE, quote=FALSE)
  write.table(pheout, pheoutfile, row.names=FALSE, quote=FALSE)
  write.table(mapout, mapoutfile, row.names=FALSE, quote=FALSE)
  
  message("Files have been converted to dlmap input format and are in working directory.")
}

