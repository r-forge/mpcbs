\name{scatter.ci}
\alias{scatter.ci}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ ~~function to do ... ~~ }
\description{
  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
scatter.ci(x, y, x.ci, y.ci, line.col = "darkgray", col = "black", ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{ ~~Describe \code{x} here~~ }
  \item{y}{ ~~Describe \code{y} here~~ }
  \item{x.ci}{ ~~Describe \code{x.ci} here~~ }
  \item{y.ci}{ ~~Describe \code{y.ci} here~~ }
  \item{line.col}{ ~~Describe \code{line.col} here~~ }
  \item{col}{ ~~Describe \code{col} here~~ }
  \item{\dots}{ ~~Describe \code{\dots} here~~ }
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
function(x,y,x.ci, y.ci, line.col="darkgray", col="black", ...){
    plot(x, y, col=col,...)
    for(i in 1:length(x)){
        segments(x[i]-x.ci[i], y[i], x[i]+x.ci[i], y[i], col=line.col)
        segments(x[i], y[i]-y.ci[i],  x[i], y[i]+y.ci[i], col=line.col)
    }
    points(x, y, col=col,...)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
