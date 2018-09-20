\name{EMMIX}
\alias{EMMIX}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Mixture Models
}
\description{
Fits finite mixture models using the Expectation-Maximization algorithm.
The mixture models can be from the (univariate and multivariate) 
normal distribution, t-distribution, skew normal distribution, and skew
t-distribution.
}
\usage{
EMMIX(dat, g, distr = "mvn", ncov = 3, clust = NULL, init = NULL,
      itmax = 1000, epsilon = 1e-06, nkmeans = 0, nrandom = 10,
      nhclust = FALSE, debug = FALSE, initloop = 20)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dat}{
  Data in the form of a \eqn{n} by \eqn{p} numeric matrix,
  where \eqn{n} is the number of observations and
  \eqn{p} the dimension of data.
}
  \item{g}{
 Number of components.
}
  \item{distr}{
String indicating the type of distribution to be fitted,
the default value is \code{"mvn"}, the normal distribution.
See Details.
}
  \item{ncov}{
A small integer indicating the type of covariance structure;
the default value is 3. See Details.
}
  \item{clust}{
A vector of integers specifying the initial partitions of the data;
the default is NULL.
}
  \item{init}{
A list containing the initial parameters for the mixture model.
See details. The default value is NULL.
}
  \item{itmax}{
An integer specifying the maximum number of iterations to apply;
the default value is 1000.
}
  \item{epsilon}{
A small number used to stop the EM algorithm loop when the relative
difference between log-likelihood at each iteration become
sufficient small; the default value is 1e-6.
}
  \item{nkmeans}{
An integer to specify the number of \code{kmeans} partitions to be used to
find the best initial values; the default value is 0
}
  \item{nrandom}{
An integer to specify the number of random partitions to be used to
find the best initial values; the default value is 10.
}
  \item{nhclust}{
A logical value to specify whether or not to use hierarchical cluster methods;
the default is FALSE. If TRUE, the Complete Linkage method will be used.
}
  \item{debug}{
A logical value, if it is TRUE, the output will be printed out; FALSE silent;
the default value is TRUE.
}
  \item{initloop}{
A integer specifying the number of initial loops when searching the best
intial partitions.
}
}
\details{
The \code{distr} determines type of the mixture distribution
which can be of 
\code{"mvn"} for a normal distribution,
\code{"mvt"} for a t-distribution, 
\code{"msn"} for a skew normal distribution, or  
\code{"mst"} for a skew t-distribution.

The \code{ncov} parameter specifies the structure of the covariance / scale 
matrix, and if it is common to components. 
Following are the available options: 
\code{ncov = 1} for a common covariance matrix,
\code{ncov = 2} for a common diagonal covariance matrix, 
\code{ncov = 3} for a general covariance/scale matrix, 
\code{ncov = 4} for a diagonal covariance/scale matrix, and
\code{ncov = 5} for \eqn{\sigma(h) * I(p)} (diagonal covariance with same
identical diagonal element values).

The \code{init} parameter can  (optionally) be used specify initial parameter
values for the EM algorithm. It should be a list containing 
following elements: 
\code{pro}, a numeric vector of the mixing proportion of each component; 
\code{mu}, a \eqn{p} by \eqn{g} matrix with each column as its 
corresponding mean; 
\code{sigma}, a three dimensional \eqn{p} by \eqn{p} by \eqn{g} array with 
its \eqn{j}th component matrix \eqn{(p, p, j)} as the covariance matrix 
for \eqn{j}th component of mixture models; 
\eqn{dof}, a vector of degrees of freedom for each component; 
\eqn{delta}, a \eqn{p} by \eqn{g} matrix with its columns
corresponding to skew parameter vectors.

When the parameter list init is given, the program ignores both initial
partition clust and automatic partition methods such as \code{nkmeans};
only when both init and clust are not available, the program uses automatic
approaches such as k-means partition method to find the best inital
values. All three automatic approaches are used to find the best initial
partition and initial values if required.

The return values include all potential parameters pro, mu, sigma, dof,
and delta, but user should not use or interpret irrelevant information
arbitrarily. For example, dof and delta for Normal mixture models.
}

\value{
\item{error}{Error code, 0 = normal exit; 1 = did not converge within
itmax iterations; 2 = failed to get the initial values; 3 = singularity.}
\item{aic}{Akaike Information Criterion (AIC).}
\item{bic}{Bayes Information Criterion (BIC).}
\item{ICL}{Integrated Completed Likelihood Criterion (ICL).}
\item{pro}{A vector of mixing proportions.}
\item{mu}{A numeric matrix with each column corresponding to the mean.}
\item{sigma}{An array of dimension \code{(p ,p, g)} with first two dimension
corresponding covari- ance matrix of each component.}
\item{dof}{A vector of degrees of freedom for each component, see Details.}
\item{delta}{A \eqn{p} by \eqn{g} matrix with each column corresponding
to a skew parameter vector.}
\item{clust}{A vector of final partition.}
\item{loglik}{The log likelihood at convergence.}
\item{lk}{A vector of log likelihood at each EM iteration.}
\item{tau}{An n by g matrix of posterior probability for each data point.}
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
 fit <- EMMIX(iris[, -5], g= 3, distr = "mvn", ncov = 3)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
