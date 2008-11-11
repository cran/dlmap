\name{dlmap.convert.cross}
\alias{dlmap.convert.cross}
\title{ Convert input files to dlmap input format }
\description{
Reads in files from all formats supported by \code{\link[qtl]{read.cross}} and converts them to dlmap input format. Creates dlmap input files in working directory.}
\usage{
dlmap.convert.cross(format = c("csv", "csvr", "csvs", "csvsr", "mm", "qtx", 
"qtlcart", "gary", "karl", "rqtl"), obj, pheobj, idname="ID", dir = "", 
genoutfile = "dlgenin.dat", pheoutfile = "dlphein.dat", 
mapoutfile = "dlmapin.dat", file, genfile, mapfile, phefile, 
chridfile, mnamesfile, pnamesfile, na.strings = c("-", "NA"),
genotypes = c("A", "H", "B", "D", "C"), alleles = c("A", "B"), ...)
}
\arguments{
 \item{format}{ See documentation for \code{\link[qtl]{read.cross}}. Also supports the 
input of an object of class \code{cross} ("rqtl" format)}
  \item{obj}{ if format="rqtl", object of class \code{cross}}
  \item{pheobj}{ Data frame or matrix containing supplementary 
spatial/environmental data}
  \item{idname}{ Unique identifier variable name; will be used to match 
phenotypic and marker data }
  \item{dir}{ Directory where input files are located; default is working 
directory}
  \item{genoutfile}{ File with genotype data output in dlmap input format }
  \item{pheoutfile}{ File with phenotype data output in dlmap input format }
  \item{mapoutfile}{ File with marker position information output in dlmap input format }
  \item{file}{ see \code{\link[qtl]{read.cross}}}
  \item{genfile}{ see \code{\link[qtl]{read.cross}}}
  \item{mapfile}{ see \code{\link[qtl]{read.cross}}}
  \item{phefile}{ see \code{\link[qtl]{read.cross}}}
  \item{chridfile}{ see \code{\link[qtl]{read.cross}}}
  \item{mnamesfile}{ see \code{\link[qtl]{read.cross}}}
  \item{pnamesfile}{ see \code{\link[qtl]{read.cross}}}
  \item{na.strings}{ see \code{\link[qtl]{read.cross}}}
  \item{genotypes}{ see \code{\link[qtl]{read.cross}}}
  \item{alleles}{ see \code{\link[qtl]{read.cross}}}
  \item{...}{ additional arguments to \code{\link[qtl]{read.cross}}}
}
\details{
Use of this function will allow automatic input into \code{\link{dlmap.asreml}} or
\code{\link{dlmap.lme}}. Otherwise input files will need to be constructed by hand. 
See below for the format if constructing input manually.

If a single set of trait values is available for each genotype, then 
phenotypic data will be input through the arguments \code{obj}, \code{file},
or \code{phefile} (depending on the file format). 
The argument \code{envobj} allows for input of phenotypic data on replicates
or additional individuals which are not necessarily genotyped. 
} 
\value{
Nothing returned. Three files output into names specified by 
\code{genoutfile}, \code{pheoutfile} and \code{mapoutfile} 
(or default values of "dlgenin.dat", "dlphein.dat", "dlmapin.dat" if 
these arguments are missing). \cr

Notation: \tabular{lll}{
\tab \dfn{n.gen} \tab is the number of genotyped individuals \cr
\tab \dfn{n.ind} \tab is the number of phenotyped individuals \cr
\tab \dfn{n.obs} \tab is the number of phenotypic observations (for a single trait) \cr
\tab \dfn{M} \tab is the number of markers \cr 
\tab \dfn{P} \tab is the number of phenotypes \cr 
\tab \dfn{C} \tab is the number of chromosomes \cr
}

Note that in general \emph{n.gen <= n.ind <= n.obs} since there may be multiple 
observations per individual, and more individuals may be phenotyped than
genotyped. Individuals which are genotyped but not phenotyped will not be
considered in analysis. \cr
\cr

\bold{Genotypes}\cr
The genotype data file \code{genoutfile} has a header row with the name of the
unique identifier variable and the marker names. This is followed by 
\emph{n.gen} rows containing the identifier and genotypes for each individual. \cr

\tabular{cccc}{
ID \tab mrk1 \tab \dots \tab mrkM \cr
1 \tab 1 \tab \dots \tab 1 \cr
\tab \dots \tab \dots \tab \cr
n.gen \tab 0 \tab \dots \tab 1 \cr
}

\bold{Phenotypes}\cr
The phenotype data file \code{pheoutfile} has a header row with the name of the
unique identifier variable and the phenotype names. This is followed by 
\emph{n.obs} rows containing the identifier and phenotypes for each sample. 
If \code{envobj} is not input, the number of rows in this file will correspond
to the number of rows in \code{genoutfile}. \cr

\tabular{cccc}{
ID \tab phename1 \tab \dots \tab phenameP \cr
1 \tab 3.7 \tab \dots \tab 2.4 \cr
\tab \dots \tab \dots \tab \cr
n.ind \tab 5.1 \tab \dots \tab 8.2 \cr
}

\bold{Linkage Map}\cr
The map data file \code{mapoutfile} has a header row with the labels shown
below for the marker ID, chromosome and position in cM. The third column can 
be omitted, in which case the marker positions will be estimated from the data. 
Note: the marker names must be in the same order as in the genotype data file. \cr

\tabular{ccc}{
MrkID \tab Chr \tab Pos \cr
mrk1 \tab 1 \tab 0.00 \cr
\tab \dots \tab \cr
mrkM \tab C \tab 100.0 \cr
}
}
\references{ Huang, BE and George, AW. Look before you leap: A new approach to QTL mapping. \emph{Manuscript in preparation}}
\author{ Emma Huang and Andrew George }
\seealso{\code{\link[qtl]{read.cross}}}
\examples{
# load dataset
data(BSdat)
data(BSphen)

\dontrun{
# convert data to dlmap format
dlmap.convert.cross(format="rqtl", obj=BSdat)

# convert data with separate phenotypic trait file
dlmap.convert.cross(format="rqtl", obj=BSdat, envobj=BSphen, idname="ID")
}
}
