\name{dlmap-internal}
\alias{dlmap.detect.asreml}
\alias{dlmap.detect.lme}
\alias{dlmap.localize.asreml}
\alias{dlmap.localize.lme}
\alias{dlmap.map.det.asreml}
\alias{dlmap.map.det.lme}
\alias{dlmap.map.loc.asreml}
\alias{dlmap.map.loc.lme}
\alias{dlmap.test.asreml}
\alias{dlmap.test.lme}
\alias{dlmap.merge}
\alias{dlmap.read}
\alias{dlmap.calc.pos}
\title{ Internal dlmap functions}
\description{
  Internal dlmap functions. These are generally not to be called by the user. 
}
\usage{
dlmap.detect.asreml(input, filestem, ...)
dlmap.detect.lme(input, filestem)

dlmap.localize.asreml(input, QTLperChr, ...)
dlmap.localize.lme(input, QTLperChr)

dlmap.map.det.asreml(input, s.chr, chrSet, prevLoc = NULL, ...)
dlmap.map.det.lme(input, s.chr, chrSet, prevLoc = NULL)

dlmap.map.loc.asreml(input, s.chr, chrSet, prevLoc = NULL, ...)
dlmap.map.loc.lme(input, s.chr, chrSet, prevLoc = NULL)

dlmap.test.asreml(input, chrSet, prevLoc = NULL, ...)
dlmap.test.lme(input, chrSet, prevLoc = NULL)

dlmap.merge(genCross, envData=NULL, step, fixpos, nperm, alpha, pindex=1, phename, idname, estmap=TRUE)

dlmap.read(genfile, phefile, mapfile, phename, filestem)

dlmap.calc.pos(exp.dat, step, fixpos)

}
\details{
The detect, localize, test and mapping functions form the basis of the 
DLMapping algorithm. There are versions for each version of dlmap (making 
use of either \command{lme} or \command{asreml} to fit the mixed models).\cr 

The merge function is used to combine the genotype and environmental data,
and also imputes missing values in the genotype data according to the Viterbi
algorithm. See \code{\link[qtl]{fill.geno}} for more details. Note that 
individuals with no phenotypic response data are omitted in the merged 
dataset. Individuals with phenotypic but no genotypic data are retained, but
the genotypes are not imputed. 
}
\references{Huang, B.E. and George, A.W. Look before you leap: A New Approach to QTL Mapping (in preparation)}
\author{ Emma Huang and Andrew George}
\keyword{internal}
