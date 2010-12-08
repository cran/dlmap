\name{dlmap-internal}
\alias{detect.asreml}
\alias{detect.lme}
\alias{localize.asreml}
\alias{localize.lme}
\alias{map.det.asreml}
\alias{map.det.lme}
\alias{map.loc.asreml}
\alias{map.loc.lme}
\alias{test.asreml}
\alias{test.lme}
\alias{calc.pos}
\alias{calc.genoprob2}
\alias{wald.test.asreml}
\alias{cintern}
\alias{dlcross.cross}
\alias{dlcross.dlmap}
\alias{dlcross.other}
\title{ Internal dlmap functions}
\description{
  Internal dlmap functions. These are generally not to be called by the user. 
}
\usage{

dlcross.cross(genobj, pheobj, idname, step, fixpos, estmap)
dlcross.dlmap(genobj, pheobj, mapobj, idname, genfile, mapfile, phefile, type, step, fixpos, estmap, ...)
dlcross.other(genobj, pheobj, mapobj, idname, genfile, mapfile, phefile)

detect.asreml(input, filestem, ...)
detect.lme(input, filestem)

localize.asreml(input, QTLperChr, ...)
localize.lme(input, QTLperChr)

map.det.asreml(input, s.chr, chrSet, prevLoc = NULL, ...)
map.det.lme(input, s.chr, chrSet, prevLoc = NULL)

map.loc.asreml(input, s.chr, chrSet, prevLoc = NULL, ...)
map.loc.lme(input, s.chr, chrSet, prevLoc = NULL)

test.asreml(input, chrSet, prevLoc = NULL, ...)
test.lme(input, chrSet, prevLoc = NULL)

calc.pos(cross, step, fixpos)
calc.genoprob2(cross, pos, error.prob=1e-04, map.function=c("haldane", "kosambi", "c-f", "morgan"))

wald.test.asreml(object, cc)
cintern(cc, tau, vrb, sigma2)
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
\references{Huang, B.E. and George, A.W. 2009. Look before you leap: A New Approach to QTL Mapping. TAG 119:899-911}
\author{ Emma Huang and Andrew George}
\keyword{internal}
