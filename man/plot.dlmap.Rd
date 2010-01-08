\name{plot.dlmap}
\alias{plot.dlmap}
\title{ Plot linkage map with detected QTL }
\description{
If the type of cross is not "other", plots the genetic linkage map for a selection of chromosomes. Indicates marker
locations, marker names, and detected QTL positions and associated flanking 
markers obtained from a \code{dlmap} fit. 
}
\usage{
\method{plot}{dlmap}(x, chr, max.dist, qcol="light blue", mcol="red", pcol="purple", \dots)
}
\arguments{
 \item{x}{object of class \code{dlmap}}
 \item{chr}{character string naming the subset of chromosomes to plot; or, if numeric, a vector of chromosome indices}
 \item{max.dist}{a numerical value in cM determining the distance the genetic map should be subsetted by}
 \item{qcol}{colour of intervals surrounding QTL (see \code{par} for colour options)}
 \item{mcol}{colour of QTL flanking markers (see \code{par} for colour options)}
 \item{pcol}{colour of QTL positions (see \code{par} for colour options)}
 \item{\dots}{arguments passed to "plot" to set up the plot region. Arguments may also be passed to "text" for the manipulation of the marker names} 
}
\details{
This function relies upon \code{link.map.cross}, which was written by Julian
Taylor for the \code{wgaim} package. It is built upon here by adding QTL 
regions and estimated positions to the map.

The function \code{plot.dlmap} provides a neat visual display of 
chromosomes. If no QTL are detected, only the linkage map will be plotted; 
otherwise detected QTL will be placed at their estimated positions 
and the intervals around them (and flanking markers) will be highlighted. 
If a subset of chromosomes are plotted and detected QTLs exist outside that 
subset a warning will be given that QTLs have been omitted from the display. 

The arguments \code{mcol}, \code{qcol} and \code{pcol} have been added for
personal colour highlighting the flanking markers, QTL regions and QTL 
positions respectively. The procedure may also be given the usual \code{col}
argument which will be passed on to the other markers. 

In order to ensure that all marker names are displayed 
without vertical overlap, the default value of the "cex" 
parameter passed to "text" should be manipulated. 
For large maps with many chromosomes, marker names and adjacent chromosomes
will overlap horizontally. In this case it is suggested that the user 
horizontally maximize the plotting window to remove overlap, or subset the
chromosomes displayed. 
} 
\value{
The genetic linkage map is plotted with shaded QTL regions, marked estimated
positions, and highlighted flanking markers. 
}
\references{ Huang, BE and George, AW. 2009. Look before you leap: A new approach to QTL mapping. TAG 119:899-911
}
\author{ Emma Huang and Andrew George; Julian Taylor }
\seealso{\code{\link[wgaim]{link.map.cross}}}

\examples{
\dontrun{
# load dataset
data(BSdat)

dl.in1 <- dlcross(format="rqtl", obj=BSdat, idname="ID", fixpos=1)
BSdl <- dlmap.asreml(object=dl.in1, phename="phenotype", filestem="BS")
BSdl <- dlmap.lme(object=dl.in1, phename="phenotype", filestem="BS")
plot(BSdl)
}
}
