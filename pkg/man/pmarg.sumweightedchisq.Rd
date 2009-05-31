\name{pmarg.sumweightedchisq}
\alias{pmarg.sumweightedchisq}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ ~~function to do ... ~~ }
\description{
  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
pmarg.sumweightedchisq(b, ALPHA, N)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{b}{ ~~Describe \code{b} here~~ }
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
function(b,ALPHA,N){
    g<-function(u){
         u^2*exp(u^2/2)/(ALPHA+exp(u^2/2))
    }
    g.moments=computeMoments(g,0);
    gnormed<-function(u){(g(u)-g.moments$psidot)/sqrt(g.moments$psidotdot)}
    theta0=sqrt(g.moments$psidotdot)/2   # need theta/sqrt(psidotdot) < 1/2 for psi(theta)<infty.
    b0=b/sqrt(N)
    THRESH=0.0001
    marg = computeTiltDirect(b0,gnormed,THRESH,theta0)
    thetabN = marg$theta*sqrt(N)
    pmarg = exp(-thetabN*b + N*log(marg$psi))/sqrt(2*pi*marg$psidotdot)
    pmarg
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
