`summary.dlcross` <- 
function(object, ...)
{
  cat("This is an object of type dlcross.\n Summary of genetic and phenotypic data:\n\n")

  if(attr(object, "type")!="other") {
	print(object$genCross)
  } else {
	gen <- as.matrix(object$dfMerged[, -c(1:object$nphe)])
     	cat(" 	Association mapping population\n\n")
	cat(" 	No. individuals: 	", ngen(object),"\n\n")
	cat("	No. phenotypes: 	", object$nphe,"\n")
	cat("	Percent phenotyped: 	", apply(object$dfMerged[,1:object$nphe], 2, function(x) return(sum(is.na(x))/length(x)*100)),"\n\n")
	cat("	No. chromosomes: 	", length(object$map), "\n\n")
	cat(" 	Total markers: 		", sum(nmrk(object)),"\n")
	cat("	No. markers:		", nmrk(object),"\n")
	cat(" 	Percent genotyped: 	", 100*sum(is.na(gen))/length(as.vector(gen)),"\n")
  }
  cat("\nThere are", ngen(object),"unique genotypes and", nphen(object), "unique phenotypes in the data.\n\n")
}
