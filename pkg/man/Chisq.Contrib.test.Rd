\name{Chisq.Contrib.test}
\alias{Chisq.Contrib.test}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ ~~function to do ... ~~ }
\description{
  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
Chisq.Contrib.test()
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
    y <- matrix(runif(100*20), ncol=20)
    T <- 5

    chpts <- matrix(runif(2*20), ncol=2)
    chpts[,2] <- rowSums(chpts)
    chpts[,1] <- as.integer(chpts[,1] * 15) + 1
    chpts[,2] <- as.integer(chpts[,2] * 15) + 3

    chpts[chpts > ncol(y)] <- ncol(y)
    
    
    res1 <- Chisq.Contrib.R(y, T, chpts)
    res2 <- Chisq.Contrib.C(y, T, chpts)

    print( range(abs(res1-res2))  )

    par(mfrow=c(1,3))
    image(res1)
    image(res2)
    image(res1-res2)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
