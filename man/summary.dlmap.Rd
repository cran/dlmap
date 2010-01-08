\name{summary.dlmap}
\alias{summary.dlmap}
\title{ Summarize dlmap object, including detected QTL }
\description{
Output summary of dlmap object, including summary of input dlcross object and
detected QTL. 
}
\usage{
\method{summary}{dlmap}(object, \dots)
}
\arguments{
 \item{object}{object of class \code{dlmap}}
 \item{\dots}{additional arguments passed to "summary" } 
}
\value{
Prints summary of input \code{dlcross} object, as well as the summary table
computed from \code{dlmap.asreml} or \code{dlmap.lme}. This includes the 
chromosome QTL are detected on, estimated positions, flanking markers, QTL
effects and standard deviations, Z-ratio and p-value. 
}
\references{ Huang, BE and George, AW. 2009. Look before you leap: A new approach to QTL mapping. TAG 119:899-911
}
\author{ Emma Huang and Andrew George; Julian Taylor }
\seealso{\code{\link[dlmap]{summary.dlcross}}, \code{\link[dlmap]{dlmap.asreml}}}

\examples{
\dontrun{
# load dataset
data(BSdat)

dl.in1 <- dlcross(format="rqtl", obj=BSdat, idname="ID", fixpos=1)
BSdl <- dlmap.asreml(object=dl.in1, phename="phenotype", filestem="BS")
BSdl <- dlmap.lme(object=dl.in1, phename="phenotype", filestem="BS")
summary(BSdl)
}
}
