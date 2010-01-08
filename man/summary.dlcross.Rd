\name{summary.dlcross}
\alias{summary.dlcross}
\title{ Summarize dlcross object }
\description{
Output summary of data stored in dlcross object, including number of genotypes, number of phenotypes, number of phenotypic variables, number of chromosomes 
and markers per chromosome.
}
\usage{
\method{summary}{dlcross}(object, \dots)
}
\arguments{
 \item{object}{object of class \code{dlcross}}
 \item{\dots}{additional arguments passed to "summary"} 
}
\value{
Prints summary of data. For experimental crosses, based off of \code{cross} 
object summary; similar output for association mapping populations. 
}
\references{ Huang, BE and George, AW. 2009. Look before you leap: A new approach to QTL mapping. TAG 119:899-911
}
\author{ Emma Huang and Andrew George; Julian Taylor }
\seealso{\code{\link[dlmap]{summary.dlmap}}, \code{\link[qtl]{summary.cross}}}

\examples{
# load dataset
data(BSdat)

dl.in1 <- dlcross(format="rqtl", genobj=BSdat, idname="ID", fixpos=1)
summary(dl.in1)
}
