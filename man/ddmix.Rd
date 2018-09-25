\name{ddmix}
\alias{ddmix}
\title{Density Functions of Mixture Models}
\description{
Calculate the density of multivariate mixture models at data points for 
each component.
}
\usage{
ddmix(dat, g, distr, mu, sigma, dof = NULL, delta = NULL, ...)
}
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
\value{
\eqn{n} by \eqn{g} matrix of log density at each data point for 
each component.
}
\references{
McLachlan GJ and Krishnan T (2008). The EM Algorithm and Extensions (2nd). 
New Jersay: Wiley.

McLachlan GJ and Peel D (2000). Finite Mixture Models. New York: Wiley.
}
\seealso{
\code{\link{ddmvn}}, \code{\link{ddmvt}}, \code{\link{ddmsn}},
\code{\link{ddmst}}. 
}
\examples{
dat <- iris[, -5]
train_ind <-  sample(1 : nrow(dat), nrow(dat) / 2, replace = FALSE)

fit <- EMMIX(dat[train_ind, ], g = 3, distr = "mvn", ncov = 3,
                   clust = iris[, 5], itmax = 100, epsilon = 0.00001)

f_dat <- ddmix(dat[train_ind, ], g = 3, distr = "mvn", ncov = 3, mu = fit$mu,
                 sigma = fit$sigma)
}
\keyword{cluster}
