\name{qchibar}
\alias{qchibar}
\title{Macro to compute quantile for chisquare mixture distribution }
\usage{
qchibar(p)
}
\arguments{
  \item{p}{ probability }
}
\description{
Reverses pvfx; i.e. given a percentile of a chisquare mixture composed of .5 chisq(df=0) and .5 chisq(df=1), computes the corresponding quantile. 
}

\author{ Emma Huang <Emma.Huang@csiro.au>}
\seealso{ \code{\link[stats]{pchisq}}, \code{\link[dlmap]{pvfx}}}
