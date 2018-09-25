\name{error.rate}
\alias{error.rate}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
 Minimum Rate of Misallocations
}
\description{
Given two vectors each corresponding to a set of categories,
this function finds the minimum number of misallocations
by rotating the categories.
}
\usage{
error.rate(clust1, clust2)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{clust1}{
A numeric or character vector of labels.
}
  \item{clust2}{
A numeric or character vector of labels same length as \code{clust1}.
}
}
\details{
Returns the minimum number of misallocations.
}
\value{
  Integer specifying the minimum number of misallocations.
}
\seealso{
\code{\link{rand.index}}
}
\examples{
clu1 <- c(1, 2, 3, 1, 1, 2, 2, 3, 3)
clu2 <- c(2, 2, 2, 1, 1, 1, 3, 3, 3)
error.rate(clu1, clu2)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
