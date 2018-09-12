\name{EMMIX2}
\alias{EMMIX2}
\title{
Mixture Models (Restricted)
}
\description{
Fits the data into the normal multivariate mixture models via the EM Algorithm, with constraints on either 1) to fix the mean values for components and then estimate the variances within components; or 2) to fix mean values for components, and specify equal variances for some components, but allow unequal variances for others.
}
\usage{
EMMIX2(dat, g, init, nvcov = 0, neq = 0, itmax = 1000, 
      epsilon = 1e-06,  debug = FALSE)
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
  \item{init}{
A list containing the initial parameters for the mixture model. 
See Details.}
  \item{nvcov}{
A vector of integer indicating which components have equal covariance 
matrices; the default value is 0, which means none of them. See Details.
}
  \item{neq}{
An integer of how many components are of the equal covariance matrices; 
the default value is 0, which means null of them. See Details.
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
  \item{debug}{
A logical value, if it is TRUE, the output will be printed out; FALSE silent; 
the default value is TRUE.
}
}
\details{
The means of the mixture have to be predefined based on prior information 
somehow and assigned to the mu of the list \code{init}. The format must be 
\eqn{p} by \eqn{g} matrix with each column corresponding to the mean of 
that component. 

\code{nvcov} has to be matched with neq, for example, if we want 
the first and third components have equal covariance matrix, 
then just specify \code{nvcov = c(1, 3)} and \code{neq = 2}.

The parameter \code{init} requires following elements: \code{pro}, 
a numeric vector of the mixing proportion of each component; 
\code{mu}|, a \eqn{p} by \eqn{g} matrix with 
each column as its corresponding mean; \code{sigma}, a three dimensional 
\code{p} by \code{p} by \code{g} array with its component matrix 
\code{sigma[„ j]} as the covariance matrix for jth component of 
mixture models.

}
\value{
\item{error} Error code, 0 = normal exit; 1 = did not converge within 
itmax iterations; 2 = failed to get the initial values; 3 = singularity.
\item{aic} Akaike Information Criterion (AIC)
\item{bic} Bayes Information Criterion (BIC)
\item{pro} A vector of mixing proportions.
\item{mu} A numeric matrix with each column corresponding to the mean.
\item{sigma} An array of dimension \code{(p ,p, g)} with first two dimension 
corresponding covari- ance matrix of each component.
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
#define the dimension of dataset
n1 <- 300
n2 <- 300
n3 <- 400;
nn <- c(n1, n2, n3)
p  <- 2
ng <- 3

#define the parameters sigma<-array(0,c(2,2,3))
for(h in 2 : 3) 
  sigma[,, h] <-diag(2) 

sigma[,, 1] <- cbind( c(1, 0.2), c(0.2, 1))
mu  <- cbind(c(4, -4), c(3.5, 4), c( 0, 0))
pro <- c(0.3, 0.3, 0.4)

distr="mvn"
# generate a data set
set.seed(111) #random seed is reset dat <- rdemmix(nn,p,ng,distr,mu,sigma)
init<-list()

init$pro<-pro 
init$mu<-mu 
init$sigma<-sigma

#for restriction 1, we have
ret1 <- EMMIX2(dat, ng, init, nvcov = 0, neq = 0, itmax = 1000, 
        epsilon = 1e-6, debug = TRUE)

#for restriction 2, we have
ret2 <- EMMIX2(dat, ng, init, nvcov = c(1, 2), neq = 2, itmax = 1000, 
        epsilon = 1e-6, debug = TRUE)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
