\name{EMMIX2}
\alias{EMMIX2}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
EMMIX2(dat, g, init, nvcov = 0, neq = 0, itmax = 1000, epsilon = 1e-06, debug = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dat}{
%%     ~~Describe \code{dat} here~~
}
  \item{g}{
%%     ~~Describe \code{g} here~~
}
  \item{init}{
%%     ~~Describe \code{init} here~~
}
  \item{nvcov}{
%%     ~~Describe \code{nvcov} here~~
}
  \item{neq}{
%%     ~~Describe \code{neq} here~~
}
  \item{itmax}{
%%     ~~Describe \code{itmax} here~~
}
  \item{epsilon}{
%%     ~~Describe \code{epsilon} here~~
}
  \item{debug}{
%%     ~~Describe \code{debug} here~~
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
function (dat, g, init, nvcov = 0, neq = 0, itmax = 1000, epsilon = 1e-06, 
    debug = FALSE) 
{
    dat <- as.matrix(dat)
    obj <- emmixfit3(dat, g, init, nvcov, neq, itmax, epsilon)
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
        ret <- obj
        if (debug) {
            msg <- paste("Restricted ", g, "- Component Multivariate Normal Mixture Model")
            cat("\n-----------------------\n\n")
            cat(msg, "\n")
            cat("\n-----------------------\n\n")
            print(obj[1:8])
        }
    }
    ret
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
