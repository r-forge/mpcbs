\name{scatter.cpfit}
\alias{scatter.cpfit}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ ~~function to do ... ~~ }
\description{
  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
scatter.cpfit(x, w, x.hat, r, alpha, platform.names)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{ ~~Describe \code{x} here~~ }
  \item{w}{ ~~Describe \code{w} here~~ }
  \item{x.hat}{ ~~Describe \code{x.hat} here~~ }
  \item{r}{ ~~Describe \code{r} here~~ }
  \item{alpha}{ ~~Describe \code{alpha} here~~ }
  \item{platform.names}{ ~~Describe \code{platform.names} here~~ }
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
function(x,w,x.hat,r,alpha, platform.names){
  K=ncol(x)
  par(mfrow=c(K-1,K-1))
        for(i in 1:(K-1)){
            for(j in (i+1):K){
               scatter.ci(x[,i], x[,j], 2*sd[,i], 2*sd[,j], cex=1.5, xlab=platform.names[[i]], ylab=platform.names[[j]], cex.lab=1.5)
               abline(alpha[j]-alpha[i]*r[j]/r[i], r[j]/r[i])
               
               points(x.hat[,i], x.hat[,j], col="red")
             #  for(a in 1:n){
             #   segments(x[a,i],x[a,j], x.hat[a,i], x.hat[a,j], col="red")
             #  }
            }
        }    
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
