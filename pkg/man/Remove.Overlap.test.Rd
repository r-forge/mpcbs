\name{Remove.Overlap.test}
\alias{Remove.Overlap.test}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ ~~function to do ... ~~ }
\description{
  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
Remove.Overlap.test()
}
%- maybe also 'usage' for other objects documented here.
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
function(){
    par(mfrow=c(1,3))

    y <- matrix(runif(2*20), ncol=2)
    y[,2] <- rowSums(y)

    i2a <- Remove.Overlap.R(y)
    i2b <- Remove.Overlap.C(y)

    

    print( sum(i2a != i2b) )

    i2a <- matrix(as.integer(i2a), nrow=1)
    i2b <- matrix(as.integer(i2a), nrow=1)

    image(i2a)
    image(i2b)
    image(i2a - i2b)
    
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
