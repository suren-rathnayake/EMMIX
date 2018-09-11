\name{EMMIX}
\alias{EMMIX}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Mixture Models 
}
\description{
fits the data into the specified multivariate mixture models via the EM 
Algorithm. Distributions (univariate and multivariate) available include 
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
the default value is \code{"mvn"}, the Normal distribution. 
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
A big integer specifying the maximum number of iterations to apply; 
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
The distribution type, determined by the distr parameter, which may take any 
one of the following values: \code{"mvn"} for a multivariate normal, 
\code{"mvt"} for a multivariate t-distribution, "msn" for a multivariate 
skew normal distribution and \code{"mst"} for a multivariate skew 
t-distribution.
The covariance matrix type, represented by the ncov parameter, 
may be any one of the following: \code{ncov = 1} for a common variance, 
\code{ncov = 2} for a common diagonal variance, \code{ncov = 3} for a 
general variance, \code{ncov = 4} for a diagonal variance, 
\code{ncov = 5} for \eqn{\sigma(h)*I(p)} (diagonal covariance with same 
identical diagonal element values).

The parameter init requires following elements: \code{pro}, a numeric vector
of the mixing proportion of each component; \code{mu}, a \eqn{p} by \eqn{g}
matrix with each column as its corresponding mean; \code{sigma}, a three 
dimensional \eqn{p} by \eqn{p} by \eqn{g} array with its \eqn{j}th 
component matrix \eqn{(p, p, j)} as the covariance matrix for \eqn{j}th
component of mixture models; \eqn{dof}, a vector of degrees of freedom for 
each component; \eqn{delta}, a \eqn{p} by \eqn{g} matrix with its columns 
corresponding to skew parameter vectors.

Since we treat the list of pro, mu, sigma, dof, and delta as a common 
structure of parameters for our mixture models, we need to include all 
of them in the initial parameter list init by default although in some 
cases it does not make sense, for example, dof and delta is not 
applicable to normal mixture model. But in most cases, the user only 
need give relevent paramters in the list.
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
\item{error} Error code, 0 = normal exit; 1 = did not converge within itmax iterations; 2 = failed to get the initial values; 3 = singularity.
\item{aic} Akaike Information Criterion (AIC)
\item{bic} Bayes Information Criterion (BIC)
\item{ICL} Integrated Completed Likelihood Criterion (ICL)
\item{pro} A vector of mixing proportions.
\item{mu} A numeric matrix with each column corresponding to the mean.
\item{sigma} An array of dimension \code{(p ,p, g)} with first two dimension 
corresponding covari- ance matrix of each component.
\item{dof} A vector of degrees of freedom for each component, see Details.
\item{delta} A \eqn{p} by \eqn{g} matrix with each column corresponding 
to a skew parameter vector.
\item{clust} A vector of final partition
\item{loglik} The log likelihood at convergence
\item{lk} A vector of log likelihood at each EM iteration
\item{tau}  An n by g matrix of posterior probability for each data point
}
\references{
%% ~put references to the literature/web site here ~
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
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (dat, g, distr = "mvn", ncov = 3, clust = NULL, init = NULL, 
    itmax = 1000, epsilon = 1e-06, nkmeans = 0, nrandom = 10, 
    nhclust = FALSE, debug = FALSE, initloop = 20) 
{
    dat <- as.matrix(dat)
    if (!is.null(init) | !missing(init)) 
        obj <- emmixfit2(dat, g, init, distr, ncov, itmax, epsilon)
    else {
        if (is.null(clust) | missing(clust)) {
            init <- try(init.mix(dat, g, distr, ncov, nkmeans, 
                nrandom, nhclust, initloop))
            if (!is.null(init)) 
                obj <- emmixfit2(dat, g, init, distr, ncov, itmax, 
                  epsilon)
            else {
                warning("not find initial values")
                obj <- list()
                obj$error = 20
            }
        }
        else {
            obj <- emmixfit1(dat, g, clust, distr, ncov, itmax, 
                epsilon, initloop)
        }
    }
    error <- obj$error
    ret <- NULL
    msg <- switch(tolower(error), `1` = paste("stopped at (did not converge within) ", 
        itmax, " iterations"), `2` = paste("density fails at initial steps! "), 
        `3` = paste("allocation fails at initial steps"), `12` = paste("density fails at estps! "), 
        `13` = paste("allocation fails at esteps"), `20` = paste("not find initials"))
    if (error > 1) {
        cat("\n-----------------------\n")
        warning("error code=", error, "\n", msg, "\n")
        cat("\n-----------------------\n\n")
    }
    if (error <= 1) {
        ICL <- getICL(dat, nrow(dat), ncol(dat), g, distr, ncov, 
            obj$pro, obj$mu, obj$sigma, obj$dof, obj$delta, obj$clust)
        ret <- obj
        ret$ICL <- ICL$ICL
        if (debug) {
            msg <- switch(tolower(distr), mvn = paste(g, "- Component Multivariate Normal Mixture Model"), 
                mvt = paste(g, "- Component Multivariate t      Mixture Model"), 
                msn = paste(g, "- Component Multivariate Skew Normal Mixture Model"), 
                mst = paste(g, "- Component Multivariate Skew-t Mixture Model"))
            cat("\n-----------------------\n\n")
            cat(msg, "\n")
            cat("\n-----------------------\n\n")
            switch(tolower(distr), mvn = print(obj[1:8]), mvt = print(obj[1:9]), 
                msn = print(obj[c(1:8, 10)]), mst = print(obj[1:10]))
            print(ICL)
            cat("\n-----------------------\n")
        }
    }
    ret
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
