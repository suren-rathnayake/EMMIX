\name{ddmvn}
\alias{ddmvn}
\alias{rdmvn}
\alias{ddmvt}
\alias{rdmvt}
\title{
The Multivariate Normal and t-Distributions
}
\description{
Density estimation and generating random variates from a
multivariate normal and t-distributions.
}
\usage{
ddmvn(dat, mean, cov)
rdmvn(n, mean, cov)
ddmvt(dat, mean, cov, nu)
rdmvt(n, mean, cov, nu)
}
%- maybe also 'usage' for other objects documented here.
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
A vector containing compoent mean.
}
  \item{cov}{
A matrix containing the compoent covariance in the case
of normal distributions, and scale matrix for t-distributions.
}
\item{nu}{
 The degrees of freedom which is a positive scaler.
}
}
\value{
Functions \code{ddmvn} and \code{ddmvt} gives the density values for
normal and \emph{t}-distributions, respectively, while \code{rdmvn} and
\code{rdmvt} generates the random numbers from normal and t-distributions,
respectively.
}
\references{
McLachlan G.J. and Krishnan T. (2008). The EM Algorithm and Extensions (2nd). New Jersay: Wiley.

McLachlan G.J. and Peel D. (2000). Finite Mixture Models. New York: Wiley.
}
\examples{
  n <- 100
  p <- 2
  mean	<- rep(0, p)
  cov	<- diag(p)
  nu	<- 3
  x	<- rdmvt(n, mean, cov, nu)
  den <- ddmvt(x, mean, cov, nu)
}

\keyword{ ~distribution }% 
\keyword{ ~multivariate }%
