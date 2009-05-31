\name{matrix.max}
\alias{matrix.max}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ ~~function to do ... ~~ }
\description{
  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
matrix.max(Z)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Z}{ ~~Describe \code{Z} here~~ }
}
\details{
  ~~ If necessary, more details than the description above ~~
}
\value{
  ~Describe the value returned
  If it is a LIST, use
  \item{comp1 }{Description of 'comp1'}
  \item{comp2 }{Description of 'comp2'}
  ...
}
\references{ ~put references to the literature/web site here ~ }
\author{ ~~who you are~~ }
\note{ ~~further notes~~ 

 ~Make other sections like Warning with \section{Warning }{....} ~
}
\seealso{ ~~objects to See Also as \code{\link{help}}, ~~~ }
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function(Z){
    # If Z has NA values treat as -Inf
    na.inds = which(is.na(Z),arr.ind=TRUE)
    Z[na.inds]=-Inf

    row.max.col=max.col(Z, ties.method=c("random", "first", "last"))
    max.row = which.max(Z[cbind(1:nrow(Z),row.max.col)])
    c(max.row,row.max.col[max.row])
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
