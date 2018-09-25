\name{rmix}
\alias{dmix}
\alias{ddmix}
\alias{rdemmix}
\alias{rdemmix2}
\alias{rdemmix3}
\title{Function for Generating Random Numbers and Density Estimations
with Finite Mixture Models}
\description{
Functions for generating random number generation and density 
estimation for finite mixture distributions of normal, \emph{t}-, 
skew normal, and skew \emph{t}-distributions.
}
\usage{
dmix(dat, model, ...)
rmix(n, model, ...)
ddmix(dat, g, distr, mu, sigma, dof = NULL, delta = NULL, ...)
rdemmix(nvect, distr, mu, sigma, dof = NULL, delta = NULL)
rdemmix2(n, distr, pro, mu, sigma, dof = NULL, delta = NULL)
rdemmix3(n, distr, pro, mu, sigma, dof = NULL, delta = NULL)
}
\arguments{
\item{dat}{
 Data in the form of a \eqn{n} by \eqn{p} numeric matrix,
 where \eqn{n} is the number of observations and
 \eqn{p} the dimension of data.
}
\item{model}{
 A list containing mixture model parameters.
}
\item{n}{
 Number of random deviates to generate.
}
\item{g}{
  Number of clusters.
}
\item{distr}{
  A string indicating the distribution; \code{"mvn"} for normal distribution,
  \code{"mvt"} for t-distribution, \code{"msn"} for skew normal distribution,
  and \code{"mst"} for skew t-distribution.
}
\item{pro}{
A vector containing of mixing proportion for the \emph{i}th 
component in the \emph{i}th element.
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
\item{nvect}{
A vector of the form \eqn{c(n_1, n_2, \dots, n_g)} where 
\eqn{n_i} is the number of variates to be generated from the
\emph{i}th component.
}
\item{...}{
  not used. Leaving it for backward compatibility.
}
}
\value{
The function \code{dmix} returns
\eqn{n} by \eqn{g} matrix of log density at each data point for 
each component.

The function \code{rmix} returns a numeric matrix with samples drawn 
in rows.

The function \code{ddmix} returns
\eqn{n} by \eqn{g} matrix of log density at each data point for 
each component.

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
\seealso{
\code{\link{ddmvn}}, \code{\link{ddmvt}}, \code{\link{ddmsn}},
\code{\link{ddmst}}. 
}
\examples{
dat <- iris[, -5]
fit <- EMMIX(dat, g = 3, distr = "mvn", ncov = 3, clust = iris[, 5], 
							itmax = 100, epsilon = 0.00001)

sim_dat <- rmix(n = 100, model = fit)
f_sim_dat <- dmix(sim_dat, model = fit)

f_dat_1 <- ddmix(sim_dat, g = 3, distr = "mvn", ncov = 3, mu = fit$mu,
                 sigma = fit$sigma)

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
\keyword{cluster}
