\name{bootstrap}
\alias{bootstrap}
\title{
Standard Errors for Mixture Model Parameters
}
\description{
Produces estimates of the standard errors for the estimated parameters
in the mixture model.
The standard errors may be assessed using parametric bootstrapping
or  nonparametric bootstrapping (i.e. by sampling with replacement).
%The standard error analysis and the bootstrap analysis of -2 log (Lambda).
}
\usage{
bootstrap(x, popPAR, B = 99, replace = TRUE, itmax = 1000, epsilon = 1e-05, ...)
}
\arguments{
  \item{x}{
  The data in the form of a \eqn{n} by \eqn{p} numeric matrix,
  where \eqn{n} is the number of observations and
  \eqn{p} the dimension of data.
}
%   \item{n}{
%   The number of observations
% }
%   \item{p}{
%   The dimension of data.
% }
%   \item{g}{
%     The number of components of the mixture model.
% }
%   \item{distr}{
% A three letter string indicating the type of distribution to be fit.
% See Details.
% }
%   \item{ncov}{
% A small integer indicating the type of covariance structure. See Details.
% }
  \item{popPAR}{
  A list containing mixture model parameters.
  % with components \code{pro}, a numeric vector containing mixing
  % proportions of each
  % component; mu, a p by g matrix with each column as its corresponding mean;
  % sigma, a three dimensional p by p by g array with its jth component matrix
  % (p,p,j) as the covariance matrix for jth component of mixture models; dof, a
  % vector of degrees of freedom for each component; delta, a p by g matrix with
  % its columns corresponding to skew parameter vectors
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
The distribution type, distr, is one of the following values: \code{"mvn"}
for a multivariate normal, \code{"mvt"} for a multivariate t-distribution,
\code{"msn"} for a multivariate skew normal distribution and \code{"mst"}
for a multivariate skew t-distribution.

The covariance matrix type, represented by the ncov parameter,
may be any one of the following: \code{ncov = 1} for a common variance,
\code{ncov = 2} for a common diagonal variance, \code{ncov = 3} for a general
variance, \code{ncov = 4} for a diagonal variance, \code{ncov = 5} for
\eqn{\sigma(h) * I(p)} (diagonal covariance with same identical diagonal
element values).

When replace is FALSE, parametric bootstrap is used;
otherwise replacement method is used.
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
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
