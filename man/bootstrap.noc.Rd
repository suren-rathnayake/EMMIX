\name{bootstrap.noc}
\alias{bootstrap.noc}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Mixture Analysis for a Range of Number of Components
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
bootstrap.noc(x, g1, g2, distr, ncov, B = 99, replace = TRUE, 
	itmax = 1000, epsilon = 1e-05)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
  The data in the form of a \eqn{n} by \eqn{p} numeric matrix,
  where \eqn{n} is the number of observations and
  \eqn{p} the dimension of data.
}
  \item{g1}{
  Minimum number of components to consider for the the mixture model.
}
  \item{g2}{
  Maximum number of components to consider for the the mixture model.
}
  \item{distr}{
 A three letter string indicating the type of distribution to be fit.
 See Details.
}
  \item{ncov}{
 A small integer indicating the type of covariance structure. See Details.
}
  \item{B}{
  The number of simulated data or replacements to be tried.
}
  \item{replace}{
  A logical value indicating whether replacement to be used.
}
  \item{itmax}{
  Maximum number of Expectation-Maximization steps for the model fitting.
}
  \item{epsilon}{
  A numeric value to test the convergence of the Expectation=Maximization
  algorithm.
  If the relative difference in log-likelihood between two (consecutive)
  EM steps falls below this values, the EM steps terminates.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
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

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
