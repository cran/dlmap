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
Version: \tab 1.0\cr
Date: \tab 2008-11-11\cr
License: \tab GPL 2\cr
}

Data for the analysis must be input in three files in the following format:
\itemize{
\item \code{genfile}: First row contains an identifier variable and marker names. 
Following rows contain genotype values for each individual
\item \code{phefile}: First row contains names of phenotypic and environmental traits. One of these traits must be an identifier for each individual. 
Following rows contain values of the traits for each individual. Note: This file may contain more observations than the file containing the marker data
\item \code{mapfile}: Contains marker names and chromosome groupings
}

The function \code{\link{dlmap.convert.cross}} is available to convert other file formats into the appropriate structure and its use is recommended to ensure that the data is properly input. 

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
\references{Huang, B.E. and George, A.W. Look before you leap: A new approach to QTL mapping. \emph{Manuscript in preparation}
}
\seealso{
\code{\link[qtl]{read.cross}}, \code{\link[asreml]{asreml}}, \code{\link[ibdreg]{ibdreg}}, \code{\link[nlme]{lme}}
}
