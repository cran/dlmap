\name{profileplot}
\alias{profileplot}
\alias{profileplot.dlmap}
\title{ Plot Wald profiles of interval mapping scans }
\description{
Plots the Wald statistic profile for each chromosome with detected QTL on first
interval mapping scan. Indicates marker locations, marker names, and detected 
QTL positions obtained from a \code{dlmap} fit.
}
\usage{
profileplot(object, \dots)

\method{profileplot}{dlmap}(object, chr, marker.names=TRUE, QTLpos=TRUE, pch=20, \dots)
}
\arguments{
 \item{object}{object of class \code{dlmap}}
 \item{chr}{character string naming the subset of chromosomes to plot; or, if numeric, a vector of chromosome indices}
 \item{marker.names}{logical value. if \code{TRUE} then marker names are plotted along the top of the plot. Defaults to \code{TRUE}}
 \item{QTLpos}{logical value. if \code{TRUE} then QTL positions are indicated with vertical lines. Defaults to \code{TRUE}}
 \item{pch}{Character to be used for points in plotting; default is solid circle}
 \item{\dots}{arguments passed to "plot" to set up the plot region. Arguments may also be passed to "text" for the manipulation of the marker names} 
}
\details{
The function \code{profileplot} provides a neat visual display of 
the Wald profile for chromosomes with detected QTL. If no QTL are detected, 
nothing will be plotted. Otherwise, the Wald profile will be plotted by cM
position of points on the interval mapping grid. Marker names will be 
displayed at the appropriate positions along the top of the plot. Vertical
lines will mark the position of detected QTL. 
}
\value{
The Wald profile from the first iteration of interval mapping scans is plotted
for all selected chromosomes, with marker names, marker locations, and QTL
positions displayed. 
}
\references{ Huang, BE and George, AW. 2009. Look before you leap: A new approach to QTL mapping. TAG 119:899-911
}
\author{ Emma Huang and Andrew George; Julian Taylor }
\examples{
\dontrun{
# load dataset
data(BSdat)

dl.in1 <- dlcross(format="rqtl", obj=BSdat, idname="ID", fixpos=1)
BSdl <- dlmap.asreml(object=dl.in1, phename="phenotype", filestem="BS")
BSdl <- dlmap.lme(object=dl.in1, phename="phenotype", filestem="BS")

profileplot(BSdl)
}
}
