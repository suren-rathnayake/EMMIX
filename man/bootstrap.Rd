\name{bootstrap}
\alias{bootstrap}
\title{
Standard Errors for Mixture Model Parameters
}
\description{
Produces estimates of the standard errors for the estimated parameters
in the mixture model.
The standard errors may be assessed using parametric 
or  nonparametric bootstrap(i.e. by sampling with replacement).
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
  A numeric value to test the convergence of the Expectation=Maximization
  algorithm.
  If the relative difference in log-likelihood between two (consecutive)
  EM steps falls below this values, the EM steps terminates.
}
\item{...}{
not used. Leaving it for backward compatibility.
}
}
\details{
When replace is FALSE, parametric bootstrap is used;
otherwise replacement method is used.
}
\value{
Standard errors of mdodel parameters.
%% ...
}
\references{
McLachlan G.J. and Krishnan T. (2008). The EM Algorithm and Extensions (2nd). New Jersay: Wiley.

McLachlan G.J. and Peel D. (2000). Finite Mixture Models. New York: Wiley.
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
  # Fit a mixture model to the iris data
  # and then obtain the standard errors of model parameters
  fit <- EMMIX(iris[, -5], g = 3, distr = "mvn", ncov = 3,
                   clust = iris[, 5], itmax = 100, epsilon = 0.00001)
  # do bootstrap (stadard error analysis)
  ## Not run:
  std <- bootstrap(iris[, -5], fit, B = 19, replace = TRUE)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~bootstrap }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~cluster }% __ONLY ONE__ keyword per line
