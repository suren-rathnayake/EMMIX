\name{bootstrap}
\alias{bootstrap}
\title{
Standard Errors for Mixture Model Parameters
}
\description{
Produces estimates of the standard errors for the estimated parameters
in the mixture model.
The standard errors may be assessed using parametric
or  nonparametric bootstrap (i.e. by sampling with replacement).
}
\usage{
bootstrap(x, model, B = 99, replace = TRUE, itmax = 1000,
          epsilon = 1e-05, ...)
}
\arguments{
  \item{x}{
  The data in the form of a \eqn{n} by \eqn{p} numeric matrix,
  where \eqn{n} is the number of observations and
  \eqn{p} the dimension of data.
}
  \item{model}{
  A list containing mixture model parameters.
}
  \item{B}{
  The number of simulated data or replacements to be tried.
}
  \item{replace}{
  A logical value indicating whether replacement to be used.
}
  \item{itmax}{
  Maximum number of Expectation-Maximization steps for the model fitting.
}
  \item{epsilon}{
  A numeric value to test the convergence of the Expectation-Maximization
  algorithm.
  If the relative difference in log-likelihood between two (consecutive)
  EM steps falls below this values, the EM steps terminates.
}
\item{...}{
not used.
}
}
\details{
When replace is FALSE, the parametric bootstrap is used;
otherwise nonparametric boostrap (replacement method) is used.
}
\value{
Standard errors of mdodel parameters.
}
\references{
McLachlan G.J. and Krishnan T. (2008). The EM Algorithm and Extensions (2nd). New Jersay: Wiley.

McLachlan G.J. and Peel D. (2000). Finite Mixture Models. New York: Wiley.
}
\examples{
set.seed(1729)
# Fit a mixture model to the iris data
# and then obtain the standard errors of model parameters
fit <- EMMIX(iris[, -5], g = 2, distr = "mvn", ncov = 3,
              itmax = 100, epsilon = 0.00001)
# do bootstrap (stadard error analysis)
## Not run:
std <- bootstrap(iris[, -5], fit, B = 19, replace = TRUE)

# print the stadard errors after rounding to three decimal points
lapply(std, round, 3)
# $se_pi
#   pi1   pi2
# 0.031 0.031

# $se_mu
#       [,1]  [,2]
# [1,] 0.058 0.042
# [2,] 0.032 0.056
# [3,] 0.086 0.023
# [4,] 0.044 0.010

# $se_sigma
# , , 1

#       [,1]  [,2]  [,3]  [,4]
# [1,] 0.058 0.027 0.066 0.022
# [2,] 0.027 0.021 0.029 0.012
# [3,] 0.066 0.029 0.080 0.032
# [4,] 0.022 0.012 0.032 0.018

# , , 2

#       [,1]  [,2]  [,3]  [,4]
# [1,] 0.020 0.018 0.008 0.003
# [2,] 0.018 0.029 0.007 0.004
# [3,] 0.008 0.007 0.007 0.002
# [4,] 0.003 0.004 0.002 0.003
}
\keyword{ ~bootstrap }
\keyword{ ~cluster }
