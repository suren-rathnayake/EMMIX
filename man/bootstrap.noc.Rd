\name{bootstrap.noc}
\alias{bootstrap.noc}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Mixture Analysis for a Range of Number of Components
}
\description{
A bootstrap approach to assess the null distribution (and hence the P-value) 
of the log likelihood ratio test (\eqn{-2 \log (\lambda)}) to test 
\eqn{H0: g = g_0} versus \eqn{H1 : g = g_0 + 1} in a given range of values
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
  A numeric value to test the convergence of the Expectation-Maximization
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
\eqn{H0: g = g_0} versus \eqn{H1 : g = g_0 + 1} tests,
and

\code{ret}: bootstrap replications of \eqn{-2 \log (\lambda)}.

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
set.seed(1089)
out <- bootstrap.noc(x = iris[, -5], g1 = 1, g2 = 4, distr = "mvn", ncov = 3, 
										B = 29, replace = TRUE, itmax = 1000, epsilon = 1e-5)


out
# $lrts
#      2 vs 1    3 vs 2   4 vs 3
# 1  346.6740  75.61674 35.67319
# 2  321.5898  77.03373 56.00523
# 3  349.9824 144.41642 67.17083
# 4  288.9804  52.79434 74.05716
# 5  348.7742 105.06670 75.59885
# 6  352.4667  77.62964 80.43432
# 7  389.3560  71.86699 76.50028
# 8  314.5501  65.77109 53.25731
# 9  356.3594  76.89652 73.97337
# 10 340.9114  84.97421 36.11725
# 11 330.8024 111.99251 36.67951
# 12 379.4971  63.52954 35.34002
# 13 344.6731 107.13048 48.21231
# 14 367.0238  89.56776 22.71701
# 15 412.0313  76.49743 88.55610
# 16 363.9714  88.78515 76.27090
# 17 386.0277  91.70892 75.69047
# 18 333.7646  96.23625 35.16773
# 19 332.6441  93.27646 77.81343
# 20 374.2174  85.79917 15.62620
# 21 305.9211  69.21048 42.08507
# 22 366.2247  64.09358 29.48164
# 23 337.7070  73.44188 84.98364
# 24 351.5601  91.14335 46.76467
# 25 375.2694  81.27428 73.37906
# 26 391.8864  39.59253 61.39599
# 27 369.1019 100.66063 94.15278
# 28 315.0528  58.13273 80.25712
# 29 321.3255 123.17055 61.43417

# $log_lk
# [1] -379.9146 -214.3547 -180.1855 -166.6645

# $pvalue
# [1] 0.24137931 0.20689655 0.06896552
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
