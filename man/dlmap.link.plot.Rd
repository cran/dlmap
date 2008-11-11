\name{dlmap.link.plot}
\alias{dlmap.link.plot}
\alias{link.map.cross}
\title{ Plot linkage map with detected QTL }
\description{
Plots the genetic linkage map for a selection of chromosomes. Indicates marker
locations, marker names, and detected QTL positions and associated flanking 
markers obtained from a \code{dlmap} fit. 
}
\usage{
dlmap.link.plot(output, chr, max.dist, marker.names=TRUE, qcol="light blue", mcol="red", pcol="purple", ...)

link.map.cross(object, chr, max.dist, marker.names = TRUE, horizontal = FALSE,  ...) 
}
\arguments{
 \item{output}{object of class \code{dlmap}}
 \item{object}{object of class \code{cross}}
 \item{chr}{character string naming the subset of chromosomes to plot; or, if numeric, a vector of chromosome indices}
 \item{max.dist}{a numerical value in cM determining the distance the genetic map should be subsetted by}
 \item{marker.names}{logical value. if \code{TRUE} then marker names are plotted alongside each chromosome on the left. Defaults to \code{TRUE}}
 \item{qcol}{colour of intervals surrounding QTL (see \code{par} for colour options)}
 \item{mcol}{colour of QTL flanking markers (see \code{par} for colour options)}
 \item{pcol}{colour of QTL positions (see \code{par} for colour options)}
 \item{horizontal}{logical value. If \code{TRUE} then map is plotted horizontally. Defaults to \code{FALSE}}
 \item{...}{arguments passed to "plot" to set up the plot region. Arguments may also be passed to "text" for the manipulation of the marker names} 
}
\details{
The function \code{link.map.cross} was written by Julian Taylor for the 
\code{wgaim} package and is built upon here by adding QTL regions and estimated
positions to the map. 

The function \code{dlmap.plot.map} provides a neat visual display of 
chromosomes. If no QTL are detected, only the linkage map will be plotted; 
otherwise detected QTL will be placed at their estimated positions 
and the intervals around them (and flanking markers) will be highlighted. 
If a subset of chromosomes are plotted and detected QTLs exist outside that 
subset a warning will be given that QTLs have been omitted from the display. 

The arguments \code{mcol}, \code{qcol} and \code{pcol} have been addedd for
personal colour highlighting the flanking markers, QTL regions and QTL 
positions respectively. The procedure may also be given the usual \code{col}
argument which will be passed on to the other markers. 

In order to ensure that all marker names are displayed 
without vertical overlap, the default value of the "cex" 
parameter passed to "text" should be mainpulated. 
For large maps with many chromosomes, marker names and adjacent chromosomes
will overlap horizontally. In this case it is suggested that the user 
horizontally maximize the plotting window to remove overlap, or subset the
chromosomes displayed. 
} 
\value{
The genetic linkage map is plotted with shaded QTL regions, marked estimated
positions, and highlighted flanking markers. 
}
\references{ Huang, BE and George, AW. Look before you leap: A new approach to QTL mapping. \emph{submitted to Genetics}
}
\author{ Emma Huang and Andrew George; Julian Taylor }
\seealso{\code{\link[wgaim]{link.map}}}

\examples{
# load dataset
data(BSdat)

\dontrun{
dlmap.convert.cross(format="rqtl", obj=BSdat)
BSdl <- dlmap.asreml(phename="phenotype", estmap=FALSE, filestem="BS")
dlmap.link.map(BSdl)
}
}
