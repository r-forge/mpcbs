\name{plotChpts}
\alias{plotChpts}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ ~~function to do ... ~~ }
\description{
  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
plotChpts(msscan.retval)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{msscan.retval}{ ~~Describe \code{msscan.retval} here~~ }
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
function(msscan.retval){

    plot(   c(0, ncol(msscan.retval$yhat)),
        c(0, nrow(msscan.retval$yhat)),
        xlab="sample index",
        ylab="SNP index",
        col="white",
        main="")
    
    for(j in 1:nrow(msscan.retval$yhat)){
        v1 <- msscan.retval$yhat[j,-1]
        v2 <- msscan.retval$yhat[j,-ncol(msscan.retval$yhat)]
        
        i2 <- which(v2 != v1)

        if(length(i2) == 0) next

        for(k in i2){
            lines(c(k, k+1, k+1, k,k), c(j,j,j+1,j+1,j), lty="solid")
        }
    }
        
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
