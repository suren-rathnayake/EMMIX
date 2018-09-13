\name{ddmix}
\alias{ddmix}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Density Functions of Mixture Models
}
\description{
Calculate the density of multivariate mixture models at data points for each 
component
}
\usage{
ddmix(dat, g, distr, mu, sigma, dof = NULL, delta = NULL, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dat}{
  Data in the form of a \eqn{n} by \eqn{p} numeric matrix,
  where \eqn{n} is the number of observations and
  \eqn{p} the dimension of data.
}
  \item{g}{
  Number of clusters.
}
  \item{distr}{
  A string indicating the distribution; \code{"mvn"} for normal distribution,
  \code{"mvt"} for t-distribution, \code{"msn"} for skew normal distribution,
  and \code{"mst"} for skew t-distribution.
}
  \item{mu}{
    A numeric mean matrix with each column containing component means. 
}
  \item{sigma}{
    An array of dimension \code{(p, p, g)} with first two dimensions 
    corresponding covariance matrix of each component.
}
  \item{dof}{
    A vector of degrees of freedom for each component.
}
  \item{delta}{
    A matrix with each column as skew parameter vector.
}
\item{...}{
    not used. Leaving it for backward compatibility.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
\eqn{n} by \eqn{g} matrix of log density at each data point for each component.
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
# subset data
dat <- iris[, -5]
train_ind <-  sample(1 : nrow(dat), nrow(dat) / 2, replace = FALSE)

fit <- EMMIX(dat[train_ind, ], g = 3, distr = "mvn", ncov = 3,
                   clust = iris[, 5], itmax = 100, epsilon = 0.00001)

ddmix(dat[train_ind, ], g = 3, distr = "mvn", ncov = 3, mu = fit$mu, 
      sigma = fit$sigma)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
