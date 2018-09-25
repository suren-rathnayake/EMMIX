\name{rdmst}
\alias{ddmst}
\alias{rdmst}
\title{The Multivariate Skew t-distribution}
\description{
Density and random generation for multivariate skew \emph{t}-distributions.
}
\usage{
  ddmst(dat, mean, cov, nu, del)
  rdmst(n, mean, cov, nu, del)
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
  A vector containing the mean parameter.
}
\item{cov}{
  A matrix containing the covariance matrix in the case
  of normal distributions, and scale matrix for \emph{t}-distributions.
}
\item{nu}{
 The degrees of freedom which is a positive scaler.
}
\item{del}{
 This is the skew parameter which is a numeric vector of length equal 
 to the dimension of the data.
}
}
\value{
\code{ddmst} gives the density values; \code{rdmst} generates the random numbers
}
\seealso{
\code{\link{rdemmix}},\code{\link{ddmvn}},\code{\link{ddmvt}},
\code{\link{ddmsn}},\code{\link{rdmvn}},\code{\link{rdmvt}},
\code{\link{rdmsn}}.
}
\examples{
n <- 100
p <- 2
mean <- rep(0, p)
cov <- diag(p)
nu <-  3
del<- c(0, 1)

y <- rdmst(n, mean, cov, nu, del)
den <- ddmst(y, mean, cov, nu, del)
}
\keyword{cluster}
\keyword{datasets}
