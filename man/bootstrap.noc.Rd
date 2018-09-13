\name{bootstrap.noc}
\alias{bootstrap.noc}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Mixture Analysis for a Range of Number of Components
}
\description{
A bootstrap approach to assess the null distribution (and hence the P-value) 
of the log likelihood ratio test (\eqn{-2 \log (\lambda)}) to test 
\eqn{H0: g = g0} versus \eqn{H1 : g = g0 + 1} in a given range of values
of \eqn{g}. 
}
\usage{
bootstrap.noc(x, g1, g2, distr, ncov, B = 99, replace = TRUE, 
							itmax = 1000, epsilon = 1e-05, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
  The data in the form of a \eqn{n} by \eqn{p} numeric matrix,
  where \eqn{n} is the number of observations and
  \eqn{p} the dimension of data.
}
  \item{g1}{
  Minimum number of components to consider for the the mixture model.
}
  \item{g2}{
  Maximum number of components to consider for the the mixture model.
}
  \item{distr}{
 A three letter string indicating the type of distribution to be fit.
 See Details.
}
  \item{ncov}{
 A small integer indicating the type of covariance structure. See Details.
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
A list containing 

\code{vlk}: log-likelihood values of fitted models 
with components from \code{g1} to \code{g2},  

\code{pvalue}: P-values for successive 
\eqn{H0: g = g0} versus \eqn{H1 : g = g0 + 1} tests,
and

\code{`ret`}: bootstrap replications of \eqn{-2 \log (\lambda)}.

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

out <- bootstrap.noc(x = iris[, -5], g1 = 3, g2 = 6, distr = "mvn", ncov = 3, 
										B = 19, replace = TRUE, itmax = 1000, epsilon = 1e-5)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
