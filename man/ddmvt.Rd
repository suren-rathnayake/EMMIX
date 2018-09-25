\name{rdmvt}
\alias{ddmvt}
\alias{rdmvt}
\title{The Multivariate t-Distribution}
\description{
Density and random generation for multivariate \emph{t}-distributions. 
}
\usage{
  ddmvt(dat, mean, cov, nu)
  rdmvt(n, mean, cov, nu)
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
\item{nu}{
 The degrees of freedom which is a positive scaler.
}
}
\value{
\code{ddmvt} gives the density values; \code{rdmvt} generates the random numbers 
}
\seealso{
\code{\link{rdemmix}},\code{\link{ddmvn}},\code{\link{ddmsn}},
\code{\link{ddmst}},\code{\link{rdmvn}},\code{\link{rdmsn}},
\code{\link{rdmst}}.
}
\examples{
n <- 100
p <- 2
mean <- rep(0, p)
cov <- diag(p)
nu <-  3

x <- rdmvt(n, mean, cov, nu)
den <- ddmvt(x, mean, cov, nu)
}
\keyword{cluster}
\keyword{datasets}
