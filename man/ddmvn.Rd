\name{ddmvn}
\alias{ddmvn}
\alias{rdmvn}
\alias{ddmvt}
\alias{rdmvt}
\alias{ddmsn}
\alias{rdmsn}
\alias{ddmst}
\alias{rdmst}
\title{
  The Multivariate Normal, t-, Skew Normal, and Skew t-Distributions
}
\description{
  Density estimation and generating random variates from a
  multivariate normal, t-, skew normal, and skew t-distributions.
}
\usage{
  ddmvn(dat, mean, cov)
  rdmvn(n, mean, cov)
  ddmvt(dat, mean, cov, nu)
  rdmvt(n, mean, cov, nu)
  ddmsn(dat, mean, cov, del)
  rdmsn(n, mean, cov, del)
  ddmst(dat, mean, cov, nu, del)
  rdmst(n, mean, cov, nu, del)
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
\item{del}{
 This is the skew parameter which is a numeric vector of length equal 
 to the dimension of the data.
}
}
\value{
  Functions \code{ddmvn}, \code{ddmvt}, \code{ddmsn} and \code{ddmst} 
  gives the density values for normal, \emph{t}-, skew normal and 
  skew \emph{t}-distributions, respectively, while 
  \code{rdmvn}, \code{rdmvt}, \code{rdmsn} and \code{rdmst} 
  generates the random these distributions, respectively.
}
\references{
  McLachlan GJ and Krishnan T (2008). The EM Algorithm and Extensions (2nd).
  New Jersay: Wiley.

  McLachlan GJ and Peel D (2000). Finite Mixture Models. New York: Wiley.
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
