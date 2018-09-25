\name{rdmvn}
\alias{ddmvn}
\alias{rdmvn}
\title{The Multivariate Normal Distribution}
\description{
Density and random generation for multivariate normal distributions. 
}
\usage{
  ddmvn(dat, mean, cov)
  rdmvn(n, mean, cov)
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
}
\value{
\code{ddmvn} gives the density values; \code{rdmvn} generates the random numbers 
}
\seealso{
\code{\link{rdemmix}},\code{\link{ddmvt}},\code{\link{ddmsn}},
\code{\link{ddmst}},\code{\link{rdmvt}},\code{\link{rdmsn}},
\code{\link{rdmst}}.
}
\examples{
n <- 100
p <- 2
mean <- rep(0, p)
cov  <- diag(p)
y <-   rdmvn(n, mean, cov)
den <- ddmvn(y, mean, cov)
}
\keyword{cluster}
\keyword{datasets}
