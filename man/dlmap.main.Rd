\name{dlmap.main}
\alias{dlmap.asreml}
\alias{dlmap.lme}
\title{ Perform DLMapping}
\description{
Fits the iterative algorithm for DLMapping. Reads in data, performs detection and localization stages and outputs summary of selected QTL effects. 
}
\usage{
dlmap.asreml(object, phename, baseModel, fixed = NULL, 
random = NULL, rcov = NULL, sparse = NULL, pedigree, 
seed = 1, n.perm = 0, multtest=c("holm", "bon"), 
alpha = 0.05, filestem = "dl", ...)

dlmap.lme(object, phename, fixed=NULL, seed=1, maxit=60, 
multtest=c("holm", "bon"), alpha=.05, filestem="dl")
}
\arguments{
  \item{object}{dlcross object created by create.dlcross}
  \item{phename}{ Response variable name}
  \item{fixed}{ A formula object specifying the fixed effects part of the 
base model, with the terms, separated by + operators, on the right of 
a ~ operator. There is no left side to the ~ expression. If no fixed effect 
is specified, the model defaults to ~1, i.e. intercept only.} 
  \item{random}{ A formula object specifying the random effects part of the 
base model, with the terms, separated by + operators, on the right of 
a ~ operator. See asreml documentation for more detail. }
  \item{rcov}{ A formula object specifying the error structure of the model, with the terms, separated by + operators, on the right of a ~ operator. See asreml documentation for more detail. }
  \item{sparse}{ A formula object specifying the fixed effects to be 
absorbed, with the terms, separated by + operators, on the right of 
a ~ operator. See asreml documentation for more detail. }
  \item{baseModel}{ An alternative to specifying \code{fixed}, \code{random}, 
\code{sparse}, and \code{rcov} separately. If a base model has already been 
fit in asreml-R for the phenotypic variation, this can be input directly}
  \item{pedigree}{ A pedigree object consisting of three columns. The first 
column is the individual ID, then the mother's ID and the father's ID. 
The name of the ID variable in the first column must match the \code{idname} variable }
  \item{seed}{ Random number seed. Default=1}
  \item{n.perm}{ Number of permutations used to get adjusted p-values at 
each iteration of detection. If \code{n.perm=0} (default) the Holm 
correction is used. }
  \item{multtest}{ Correction used for multiple testing. If n.perm>0 will use 
  permutation, but otherwise can choose between Holm and Bonferroni}
  \item{alpha}{ Significance level for testing }
  \item{filestem}{ Stem to add to names of any files generated in DL Mapping 
process. Default="dl"}
  \item{maxit}{ Maximum number of iterations to attempt for convergence of \code{lme}}
  \item{...}{ additional arguments to asreml}
}
\details{
There are two versions of this function, which use different engines to 
fit the linear mixed models which form the framework of the algorithm.
\code{dlmap.asreml} provides a much more general implementation of the
DLMapping algorithm and is the preferred method of analysis. 
\code{dlmap.lme} is more restricted in its capabilities, in 
that it cannot model random effects or covariance structure, cannot handle 
more than 200 markers, and only allows for a single phenotypic observation
per genotype. Also, permutation has not been implemented for this function
because it is very slow. However, \code{dlmap.lme} will fit the basic 
algorithm and is useful should a license for ASReml not be available. 

In \code{dlmap.asreml}, there are two options for specifying the model for 
phenotypic variation. 
The individual model components can either be input directly as they would be 
in an ASReml call, or a previous model (\code{baseModel}) output from ASReml 
can be input and the components will be retrieved from it. The latter 
formulation may be useful if prior phenotypic modelling has taken place. Note
that in either case, variables appearing in the rcov statement must be 
ordered appropriately in the dataset. For example, if 
\command{rcov=~ar1(Column):ar1(Row)} the data must be sorted as 
\emph{Row} within \emph{Column}. 

Missing values in \command{asreml} are replaced with zeros, so it is important
to centre the covariate in question. This is done for all genotypes within the
\command{dlmap.asreml} function. Thus individuals with phenotypic but not
genotypic data, which play important roles in field trials, may be included
safely. For \command{dlmap.lme} these individuals cannot be included, so the
default behavior is to omit observations with missing values. 

It is recommended that \code{no.perm} be set to 0 for initial exploratory 
analysis, as the permutation analysis may be lengthy. The Holm 
correction is used to adjust for the number of chromosomes under 
consideration at each detection stage. While this is a conservative 
measure it seems to perform well in practice.

Two files are output with names set by the argument \code{filestem}, which
has a default value of "dl". The 
file "filestem.trace" contains ASReml licensing and likelihood convergence
output which otherwise would be dumped to the screen and possibly obscure
other messages. Errors, warnings and other messages will still appear on the 
screen. Some warnings which appear may be passed through from an ASReml call
and output on exit. These may generally be ignored. This file is not created
if \code{dlmap.lme} is used.  

The file "filestem.det.log" is a record of iterations in the detection stage. 
For each iteration the REMLRT testing for genetic variation on each chromosome
is output, along with adjusted p-values, genomewide threshold and markers 
selected as fixed effects. The p-values are corrected for the
number of chromosomes tested either by the Holm correction or by 
permutation. If the number of permutations (\code{n.perm}) is greater than 
0, then for the Xth iteration an additional file "filestem.permX" will be 
created which contains the test statistics for the permuted datasets. 
See the accompanying vignette for an example of how to interpret the ".det.log" file. 

}
\value{
  \item{Summary}{Table with one row per QTL detected, columns for which 
chromosome the QTL is on, its position (cM), flanking markers, additive 
(dominant) effect and standard deviation, Z-ratio and p-value.}
  \item{no.qtl }{ Total number of QTL detected on all chromosomes}
  \item{final.model }{Object of class \code{asreml} for final model 
containing all terms in the base model, as well as effects for every QTL 
detected at the appropriate locations. No random effects for markers are fit}
  \item{profile}{If QTL are detected on \emph{C} chromosomes, this is a list 
with \emph{C} elements, each a matrix with 2 rows and a column for each 
position on the chromosome. The first row contains the cM position; the second 
row contains the Wald statistic for the model fit in the localization stage}
  \item{input}{Original dlcross input to the analysis}
}
\references{ Huang, BE and George, AW. 2009. Look before you leap: A new approach to QTL mapping. TAG 119:899-911}
\author{ Emma Huang and Andrew George }
\seealso{ \code{\link[dlmap]{dlcross}}, \code{\link[dlmap]{plot.dlmap}}, \code{\link[dlmap]{summary.dlmap}} }
\examples{
\dontrun{
data(BSdat)
data(BSphe)

# Convert cross object to DL Mapping format
dl.in1 <- dlcross(format="rqtl", genobj=BSdat, idname="ID", fixpos=1)

# Analyze data
BSdl <- dlmap.asreml(object=dl.in1, phename="phenotype", filestem="BS")
BSdl <- dlmap.lme(object=dl.in1, phename="phenotype", filestem="BS")

# With additional phenotypic data
dl.in2 <- dlcross(format="rqtl", genobj=BSdat, pheobj=BSphe, idname="ID", step=5)
BSph <- dlmap.asreml(object=dl.in2, phename="phenotype", env=TRUE, random=~Block)
}
}
