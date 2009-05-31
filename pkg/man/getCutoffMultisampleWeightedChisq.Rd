\name{getCutoffMultisampleWeightedChisq}
\alias{getCutoffMultisampleWeightedChisq}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ ~~function to do ... ~~ }
\description{
  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
getCutoffMultisampleWeightedChisq(pval, m, delta, win, N, alpha)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{pval}{ ~~Describe \code{pval} here~~ }
  \item{m}{ ~~Describe \code{m} here~~ }
  \item{delta}{ ~~Describe \code{delta} here~~ }
  \item{win}{ ~~Describe \code{win} here~~ }
  \item{N}{ ~~Describe \code{N} here~~ }
  \item{alpha}{ ~~Describe \code{alpha} here~~ }
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
function(pval,m,delta,win,N,alpha){
    
    cat("Computing threshold for weighted chi-square...\n")
    THRES = 0.1*pval
    currb = 1
    prevsmallerb = currb
    prevlargerb = 200
    currpval=1
    
    while( abs(currpval-pval)>THRES ){
#        cat("pval =", pval, ", currpval = ",currpval,", THRES=", THRES,".\n",sep="")
        
        if( currpval>pval){
            # need to increase b.
            prevsmallerb = currb
            currb = currb+(prevlargerb-currb)/2
        } else {
            # need to decrease b.
            prevlargerb = currb
            currb = currb - (currb-prevsmallerb)/2
        }
    
#        cat("currb = ",currb,"\n")
        currpval = pvalueMultisampleWeightedChisq(currb,m,delta,win,alpha,N);
    }    
    currb
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
