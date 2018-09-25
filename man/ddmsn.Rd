\name{rdmsn}
\alias{ddmsn}
\alias{rdmsn}
\title{The Multivariate Skew Normal Distribution}
\description{
Density and random generation for multivariate skew normal 
distributions. 
}
\usage{
  ddmsn(dat, mean, cov, del)
  rdmsn(n, mean, cov, del)
}
\arguments{
\item{dat}{
  Data in the form of a \eqn{n} by \eqn{p} numeric matrix,
  where \eqn{n} is the number of observations and
  \eqn{p} the dimension of data. 
}
\item{n}{
  Number of observations (to be drawn).
}
\item{mean}{
  A vector containing mean parameter.
}
\item{cov}{
  A matrix containing the covariance matrix in the case
  of normal distributions, and scale matrix for \emph{t}-distributions.
}
\item{del}{
 This is the skew parameter which is a numeric vector of length equal 
 to the dimension of the data.
}
}
\value{
\code{ddmsn} gives the density values; \code{rdmsn} generates the random numbers 
}
\seealso{
\code{\link{rdemmix}},\code{\link{ddmvn}},\code{\link{ddmvt}},
\code{\link{ddmst}},\code{\link{rdmvn}},\code{\link{rdmvt}},
\code{\link{rdmst}}.
}
\examples{
n <- 100
p <- 2
mean <- rep(0, p)
cov  <- diag(p)
del<- c(0, 1)

y   <- rdmsn(n, mean, cov, del)
den <- ddmsn(y, mean, cov, del)
}
\keyword{cluster}
\keyword{datasets}
