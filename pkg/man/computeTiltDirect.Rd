\name{computeTiltDirect}
\alias{computeTiltDirect}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ ~~function to do ... ~~ }
\description{
  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
computeTiltDirect(b, g, THRESH, theta0)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{b}{ ~~Describe \code{b} here~~ }
  \item{g}{ ~~Describe \code{g} here~~ }
  \item{THRESH}{ ~~Describe \code{THRESH} here~~ }
  \item{theta0}{ ~~Describe \code{theta0} here~~ }
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
function(b,g,THRESH,theta0) {

    theta = theta0 # if start theta at 0, then dh(theta)=0 at the first step for g(u)=u^2.
    prevtheta=Inf
    prevprevtheta=Inf
    thetarec = theta
    
    while( abs(theta-prevtheta)>THRESH && abs(theta-prevprevtheta)>THRESH ){
        
        g.moments<-computeMoments(g,theta)
        htheta=g.moments$psidot-b
        dhtheta = g.moments$psidotdot

        # cat("theta=", theta," htheta=",htheta,".\n",sep="")
        thetarec= c(thetarec, theta)
        prevprevtheta=prevtheta
        prevtheta=theta
        theta = prevtheta - htheta/dhtheta
    }

    theta=prevtheta
    psi = g.moments$psi

    list(theta=theta,psi=psi,psidot=g.moments$psidot,psidotdot=g.moments$psidotdot)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
