\name{pvfx}
\alias{pvfx}
\title{ Macro to compute p-value for testing variance component = 0}
\usage{
pvfx(test)
}
\arguments{
  \item{test}{ Test statistic for which to compute p-value}
}
\description{
Given a value > 0, computes 1 - (percentile of a chisquare mixture composed of .5 chisq(df=0) and .5 chisq(df=1)).
}
\author{ Emma Huang and Andrew George }
\seealso{ \code{\link[stats]{pchisq}}, \code{\link[dlmap]{qchibar}}}
