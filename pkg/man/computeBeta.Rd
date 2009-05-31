\name{computeBeta}
\alias{computeBeta}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ ~~function to do ... ~~ }
\description{
  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
computeBeta(ALPHA)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{ALPHA}{ ~~Describe \code{ALPHA} here~~ }
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
function(ALPHA){
    w<-function(u){
         exp(u^2/2)/(ALPHA+exp(u^2/2))
    }
    integrand<-function(u){
        u^2*w(u)^2*(2*u^4*w(u)*(1-w(u))+5*u^2*w(u)-3*u^2-2)*dchi(u,1)
    }
    
#    us=seq(0,10,0.1)
#    ys=integrand(us)
#    plot(us,ys)
    
    numerator=integrate(integrand,lower=0,upper=10)
    
    g<-function(u){
         u^2*exp(u^2/2)/(ALPHA+exp(u^2/2))
    }
    g.moments=computeMoments(g,0)
    numerator$value/(2*g.moments$psidotdot)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
