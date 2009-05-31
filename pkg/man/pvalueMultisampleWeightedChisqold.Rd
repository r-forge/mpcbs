\name{pvalueMultisampleWeightedChisqold}
\alias{pvalueMultisampleWeightedChisqold}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ ~~function to do ... ~~ }
\description{
  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
pvalueMultisampleWeightedChisqold(b, m, delta, win, ALPHA, N)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{b}{ ~~Describe \code{b} here~~ }
  \item{m}{ ~~Describe \code{m} here~~ }
  \item{delta}{ ~~Describe \code{delta} here~~ }
  \item{win}{ ~~Describe \code{win} here~~ }
  \item{ALPHA}{ ~~Describe \code{ALPHA} here~~ }
  \item{N}{ ~~Describe \code{N} here~~ }
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
function(b,m,delta,win,ALPHA,N){
    delta1 = min(1, win/m);
#    if(msscan.debug.trace){
#        print(paste("m = ", m, "; b = ", b, "; delta = ", delta, "; delta1 = ", delta1))
#    }

    beta=computeBeta(ALPHA)

    integrand<-function(u){
        vu(sqrt(2)*b/(sqrt(m)*sqrt(u*(1-u))))^2/(u^2*(1-u))
    }
    
    integral =  integrate(integrand,lower=delta,upper=delta1)
    pmarg = pmarg.sumweightedchisq(b,ALPHA,N)

    pval=b^3*beta^2*pmarg*integral$value
    pval
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
