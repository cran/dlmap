\name{dlcross}
\alias{dlcross}
\title{ Constructor function for dlcross format}
\description{
Reads in files from all formats supported by \code{\link[qtl]{read.cross}} and converts them to dlmap input format. Also supports other formats for objects}
\usage{
dlcross(format = c("csv", "csvr", "csvs", "csvsr", "mm", "qtx", 
"qtlcart", "gary", "karl", "rqtl", "dlmap", "other"), genobj, pheobj, 
mapobj, idname="ID", dir = "", file, genfile, mapfile, phefile, 
chridfile, mnamesfile, pnamesfile, na.strings = c("-", "NA"),
genotypes = c("A", "H", "B", "D", "C"), alleles = c("A", "B"), type, step=0,
fixpos=0, estmap=TRUE, ...)
}
\arguments{
 \item{format}{ See documentation for \code{\link[qtl]{read.cross}}. Also supports the 
input of an object of class \code{cross} ("rqtl" format), "dlmap" format 
(described below) and continuous data for association mapping populations 
("other" format)}
  \item{genobj}{ if format="rqtl", object of class \code{cross}; 
if format="dlmap" or "other" data frame containing genotypes}
  \item{pheobj}{ Data frame or matrix containing supplementary 
phenotypic or environmental data}
  \item{mapobj}{ if format="other" or "dlmap", 2/3-columned data-frame 
containing marker names, chromosomes and marker positions}
  \item{idname}{ Unique identifier variable name; will be used to match 
phenotypic and marker data }
  \item{dir}{ Directory where input files are located; default is working 
directory}
  \item{file}{ see \code{\link[qtl]{read.cross}}}
  \item{genfile}{ see \code{\link[qtl]{read.cross}}; if format="dlmap" see below}
  \item{mapfile}{ see \code{\link[qtl]{read.cross}}; if format="dlmap" see below}
  \item{phefile}{ see \code{\link[qtl]{read.cross}}; if format="dlmap" see below}
  \item{chridfile}{ see \code{\link[qtl]{read.cross}}}
  \item{mnamesfile}{ see \code{\link[qtl]{read.cross}}}
  \item{pnamesfile}{ see \code{\link[qtl]{read.cross}}}
  \item{na.strings}{ see \code{\link[qtl]{read.cross}}}
  \item{genotypes}{ see \code{\link[qtl]{read.cross}}}
  \item{alleles}{ see \code{\link[qtl]{read.cross}}}
  \item{type}{Type of experimental cross; only necessary if format="dlmap" or 
"other"; if not input, will be assumed that data is from unrelated individuals}
  \item{step}{ Step size for localization stage, i.e. if \code{step=2}, grid 
of positions spaced 2 cM apart are considered for QTL locations. If 
\code{step=0} (default) positions are only located at markers.}
  \item{fixpos}{ Alternative to specifying a step size - if \code{fixpos=2}, 
2 evenly spaced positions between each marker are considered as QTL locations. 
If \code{fixpos=0} (default) positions are only located at markers. }
  \item{estmap}{Flag for whether to re-estimate the linkage map. Cannot be done 
for format="other"}
  \item{...}{ additional arguments to \code{\link[qtl]{read.cross}}}
}
\details{
This function constructs \code{dlcross} objects for input into 
\code{\link{dlmap.asreml}} or \code{\link{dlmap.lme}}. The \code{format} 
argument allows for a wide range of input data formats, including all those
compatible with \code{\link[qtl]{read.cross}}. In addition, format="dlmap"
takes three files or objects in the following form. 

\itemize{
\item \code{genfile}: First row contains an identifier variable and marker names. 
Following rows contain genotype values for each individual
\item \code{phefile}: First row contains names of phenotypic and environmental traits. One of these traits must be an identifier for each individual. 
Following rows contain values of the traits for each individual. Note: This file may contain more observations than the file containing the marker data
\item \code{mapfile}: Contains marker names, chromosome groupings and positions
}

If a single set of trait values is available for each genotype, then 
phenotypic data will be input through the arguments \code{genobj}, \code{file},
or \code{phefile} (depending on the file format). 
The argument \code{pheobj} allows for input of phenotypic data on replicates
or additional individuals which are not necessarily genotyped. 

Choosing format="other" allows for association mapping populations to be 
analyzed, and in this case the data can be input as with format="dlmap", 
but a genetic linkage map is not required. Hence the mapfile only needs the
first two columns of marker names and chromosome groupings. 

} 
\value{ Object with class "dlcross" which can be input to dlmap.asreml or 
dlmap.lme. Contains the following elements:
\item{dfMerged}{Data frame to be used in dlMapping analysis}
\item{map}{Original genetic map}
\item{nphe}{Number of phenotypic traits}
\item{loc}{A flag for whether to run the localization stage}
\item{idname}{the ID name input to the function}
\item{mapp}{If format is not "other", genetic map augmented by imputed genotypes at grid of positions}
\item{genCross}{if format is not "other", rqtl cross object. CAUTION: if 
there are replicates, i.e., more phenotypic data than genotypic, this object
will not contain all of the phenotypic data for the sample}
}
\seealso{\code{\link[dlmap]{summary.dlcross}}, \code{\link[dlmap]{plot.dlcross}}}
\references{ Huang, BE and George, AW. 2009. Look before you leap: A new approach to QTL mapping. TAG 119:899-911}
\author{ Emma Huang and Andrew George }
\examples{
# load dataset
data(BSdat)
data(BSphe)

# convert data to dlmap format
dl.in1 <- dlcross(format="rqtl", genobj=BSdat, idname="ID", fixpos=1)

# convert data with separate phenotypic trait file
dl.in2 <- dlcross(format="rqtl", genobj=BSdat, pheobj=BSphe, idname="ID", step=5)
}
