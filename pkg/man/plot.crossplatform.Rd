\name{plot.crossplatform}
\alias{plot.crossplatform}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ ~~function to do ... ~~ }
\description{
  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
plot.crossplatform(pos, y, anchor, chpts = NULL, ranking = NULL, yhat = NULL, platform.names, num.panels = NA, xlim = NULL, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{pos}{ ~~Describe \code{pos} here~~ }
  \item{y}{ ~~Describe \code{y} here~~ }
  \item{anchor}{ ~~Describe \code{anchor} here~~ }
  \item{chpts}{ ~~Describe \code{chpts} here~~ }
  \item{ranking}{ ~~Describe \code{ranking} here~~ }
  \item{yhat}{ ~~Describe \code{yhat} here~~ }
  \item{platform.names}{ ~~Describe \code{platform.names} here~~ }
  \item{num.panels}{ ~~Describe \code{num.panels} here~~ }
  \item{xlim}{ ~~Describe \code{xlim} here~~ }
  \item{\dots}{ ~~Describe \code{\dots} here~~ }
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
function(pos,y,anchor,chpts=NULL,ranking=NULL,yhat=NULL,platform.names,num.panels=NA, xlim=NULL,...){
    K = length(y)
    
    if(is.na(num.panels)) num.panels=K
    
    par(mfrow=c(num.panels,1))
    if(is.null(xlim)){
        xlim=0
        
        xlim=c(min(anchor$merged.pos),xmax=max(anchor$merged.pos))
    } 
        
    for(k in 1:K){
        ymax = max(y[[k]]) + 0.2* max(y[[k]])
        ymin = min(y[[k]]) - 0.2*abs(min(y[[k]]))
    
        plot(pos[[k]],y[[k]], xlim=xlim,ylim=c(ymin,ymax),main=platform.names[[k]],ylab="",xlab="base position", ...)
        if(!is.null(yhat)) lines(pos[[k]],yhat[[k]], col="red",lwd=2)
        if(length(chpts)>0) {
            segments(anchor$merged.pos[chpts],rep(ymin,length(chpts)),anchor$merged.pos[chpts],rep(ymax,length(chpts)),col="black",lwd=1)
            
            if(!is.null(ranking)){
                labels = format(ranking)
                text(anchor$merged.pos[chpts], rep(max(y[[k]]),length(chpts)), labels,adj = c(0,-1))
            }
        }
    }    
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
