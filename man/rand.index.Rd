\name{rand.index}
\alias{rand.index}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Rand Index
}
\description{
Computes the Rand index and the adjusted Rand index. }
\usage{
rand.index(LabelA, LabelB)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{LabelA}{
    A numeric or character vector of labels.
}
  \item{LabelB}{
  A numeric or character vector of labels same length as \code{LabelA}.
}
}
\details{
Measures the agreement between two sets of partitions. 
The upper bound of 1 implies perfect agreement. 
The expected value is zero if the partitions are random.
}
\value{
Scaler specifying how closely two partitions agree.   
}
\references{
Hubert L, and Arabie P (1985). Comparing Partitions.
\emph{Journal of the Classification} \strong{2}, 193--218.
}
\seealso{
\code{\link{error.rate}}
}
\examples{
clu1 <- c(1, 2, 3, 1, 1, 2, 2, 3, 3)
clu2 <- c(2, 2, 2, 1, 1, 1, 3, 3, 3)
rand.index(clu1, clu2)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
