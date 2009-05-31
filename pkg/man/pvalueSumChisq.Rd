\name{pvalueSumChisq}
\alias{pvalueSumChisq}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ ~~function to do ... ~~ }
\description{
  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
pvalueSumChisq(b.normed, T, delta, T0, N)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{b.normed}{ ~~Describe \code{b.normed} here~~ }
  \item{T}{ ~~Describe \code{T} here~~ }
  \item{delta}{ ~~Describe \code{delta} here~~ }
  \item{T0}{ ~~Describe \code{T0} here~~ }
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
function(b.normed,T,delta,T0,N){

    b=sqrt(b.normed*sqrt(2*N)+N)
    
    delta1 = min(1, T0/T)
    CbN = 1-(N-1)/b^2
    
    integrand<-function(u){
        vu(b*CbN/(sqrt(T)*sqrt(u*(1-u))))^2/(u^2*(1-u))
    }
    
    integral =  integrate(integrand,lower=delta,upper=delta1)
    
    pval=b^3*dchi(b,N)*CbN^3*2^{-2}*integral$value
    pval
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
