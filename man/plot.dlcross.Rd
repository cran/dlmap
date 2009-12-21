\name{plot.dlcross}
\alias{plot.dlcross}
\title{ Plot diagnostics summarizing data}
\description{
If the type of cross is not "other", will plot the genetic map for the cross. Will also plot a histogram or barplot of the first few phenotypic variables. 
}
\usage{
\method{plot}{dlcross}(x, chr, pheno.col, \dots)
}
\arguments{
 \item{x}{object of class \code{dlcross}}
 \item{chr}{character string naming the subset of chromosomes to plot; or, if numeric, a vector of chromosome indices}
 \item{pheno.col}{vector of phenotypes for which to show barplot/histogram}
 \item{\dots}{arguments passed to "plot" to set up the plot region. Arguments may also be passed to "text" for the manipulation of the marker names} 
}
\details{
Will only plot a maximum of three phenotypic variables; pheno.col can be 
used to select which are plotted. 
} 
\references{ Huang, BE and George, AW. 2009. Look before you leap: A new approach to QTL mapping. TAG 119:899-911
}
\author{ Emma Huang and Andrew George; Julian Taylor }
\seealso{\code{\link[qtl]{plot.cross}}}

\examples{
# load dataset
data(BSdat)

\dontrun{
dl.in1 <- dlcross(format="rqtl", obj=BSdat, idname="ID", fixpos=1)
plot(dl.in1)
}
}
