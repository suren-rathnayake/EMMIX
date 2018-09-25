\name{rdemmix}
\alias{rdemmix}
\alias{rdemmix2}
\alias{rdemmix3}
\title{Simulate Data Using Mixture Models}
\description{
Generate random number from specified mixture models,
including univariate and multivariate normal distribution,
\emph{t}-distribution, skew normal distribution,
and skew \emph{t}-distribution.
}
\usage{
rdemmix(nvect, distr, mu, sigma, dof = NULL, delta = NULL)
rdemmix2(n, distr, pro, mu, sigma, dof = NULL, delta = NULL)
rdemmix3(n, distr, pro, mu, sigma, dof = NULL, delta = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{nvect}{
A vector of the form \eqn{c(n_1, n_2, \dots, n_g)} where 
\eqn{n_i} is the number of variates to be generated from the
\emph{i}th component.
}
\item{distr}{
A three letter string indicating the distribution type.
}
\item{n}{
Number of observations.
}
\item{pro}{
A vector containing of mixing proportion for the \emph{i}th 
component in the \emph{i}th element.
}
\item{mu}{
A numeric matrix of which the \emph{i}th column corresponds to the mean 
of the \emph{i}th component. 
}
\item{sigma}{
A three-dimensional array of where \code{sigma[,, i]} 
corresponds to the covariance/scale matrix of the \emph{i}th component. 
}
\item{dof}{
A vector containing the degrees of freedom paraneter of the \emph{i}th 
component in the \emph{i}th element. Defined only for
\emph{t}- and the skew \emph{t}-distributions. 
}
\item{delta}{
A numeric matrix of which the \emph{i}th column corresponds to the skew  
parameter vector of the \emph{i}th component. Defined only 
skew normal and skew \emph{t}-distributions.
}
}
\details{
The distribution type, determined by the \code{distr} parameter,
which may take any one of the following values:
"mvn" for a multivariate normal, "mvt" for a multivariate
\emph{t}-distribution, "msn" for a multivariate skew normal
distribution and "mst" for a multivariate skew \emph{t}-distribution.
}
\value{
Both \code{rdemmix} and \code{rdemmix2} return an \eqn{n} by \eqn{p} 
numeric matrix with rows corresponding to the random variates 
generated, where as \code{rdemmix3} returns a list which, 
in addition to the random variates generated contain the label,
\code{cluster}, from which each variate is generated.
}
\references{
McLachlan GJ and Krishnan T (2008). The EM Algorithm and Extensions (2nd).
New Jersay: Wiley.

McLachlan GJ and Peel D (2000). Finite Mixture Models. New York: Wiley.
}
\examples{
#specify the dimension of data, and number of clusters
#the number of observations in each cluster
n1 <- 300
n2 <- 300
n3 <- 400
nn <- c(n1, n2, n3)

p <- 2
g <- 3

#specify the distribution
distr <- "mvn"

#specify mean and covariance matrix for each component
sigma <- array(0, c(2, 2, 3))
for (h in 2 : 3)
   sigma[,,h] <- diag(2)

sigma[,, 1] <- cbind(c(1, -0.1), c(-0.1, 1))
mu  <- cbind(c(4, -4), c(3.5, 4), c(0, 0))

#reset the random seed
set.seed(111)
#generate the dataset
dat <- rdemmix(nn, distr, mu, sigma)

# alternatively one can use
pro   <- c(0.3, 0.3, 0.4)
n <- 1000
set.seed(111)
dat <- rdemmix2(n, distr, pro, mu, sigma)
plot(dat)

# and
set.seed(111)
dobj <- rdemmix3(n, distr, pro, mu, sigma)
plot(dobj$data)

#other distributions such as "mvt","msn", and "mst".
#t-distributions
dof    <- c(3, 5, 5)
dat <- rdemmix2(n, "mvt", pro, mu, sigma, dof)
plot(dat)

#Skew Normal distribution
delta <- cbind(c(3, 3), c(1, 5), c(-3, 1))
dat <- rdemmix2(n, "msn", pro, mu, sigma, delta = delta)
plot(dat)

#Skew t-distribution
dat <- rdemmix2(n, "mst", pro, mu, sigma, dof, delta)
plot(dat)
}
\keyword{ ~distribution }%
\keyword{ ~multivariate }% 
