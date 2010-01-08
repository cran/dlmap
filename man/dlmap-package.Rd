\name{dlmap-package}
\alias{dlmap-package}
\alias{dlmap}
\docType{package}
\title{
DLMapping for QTL detection
}
\description{
QTL mapping in a mixed model framework with separate detection and 
localization stages. The former detects the number of QTL on each chromosome 
based on the genetic variation due to the grouped markers on the chromosome, 
while the latter stage uses this information to determine the most likely QTL 
positions. The mixed model can accommodate general fixed and random effects, 
including spatial effects in field trials and random pedigree effects. 
}
\details{
\tabular{ll}{
Package: \tab dlmap\cr
Type: \tab Package\cr
Version: \tab 1.04\cr
Date: \tab 2010-01-08\cr
License: \tab GPL 2\cr
}

The function \code{\link{dlcross}} is a constructor function for \code{dlcross}
objects to be input to \code{\link{dlmap.asreml}} or \code{\link{dlmap.lme}}. 
It will read in files in several different formats, including all those 
compatible with \code{\link[qtl]{read.cross}} and two new formats to 
accommodate association mapping populations and designs with extensive 
phenotypic data. 

The primary function is \code{\link{dlmap.asreml}}, which performs the iterative 
algorithm to detect and position QTL on all chromosomes with significant 
genetic variation. This can accomodate sophisticated mixed models for phenotypic
variation in addition to the genetic modeling. 

Because ASReml-R is proprietary, we provide additional functions for DLMapping using the \pkg{nlme} package to fit mixed models. However, the function
\code{\link{dlmap.lme}} does not provide as many features as that based on ASReml-R. 

The vignette included in this package gives more background on the methodology, input file structure, and examples of how to use each of the important functions
in the package. 

}
\author{
Emma Huang and Andrew George \cr 
Maintainer: Emma Huang <Emma.Huang@csiro.au>
}
\references{Huang, B.E. and George, A.W. 2009. Look before you leap: A new approach to QTL mapping. TAG 119:899-911
}
\seealso{
\code{\link[ibdreg]{ibdreg}}, \code{\link[nlme]{lme}}
}
