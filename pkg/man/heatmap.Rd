\name{heatmap}
\alias{heatmap}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ ~~function to do ... ~~ }
\description{
  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
heatmap(y, chrom = NULL, loc = NULL, horizontal = TRUE, pos.string = "Probe order", ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{y}{ ~~Describe \code{y} here~~ }
  \item{chrom}{ ~~Describe \code{chrom} here~~ }
  \item{loc}{ ~~Describe \code{loc} here~~ }
  \item{horizontal}{ ~~Describe \code{horizontal} here~~ }
  \item{pos.string}{ ~~Describe \code{pos.string} here~~ }
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
function(y,chrom=NULL,loc=NULL, horizontal=TRUE, pos.string="Probe order",...){
    nsamples = ncol(y)
    if(is.null(loc)) loc=1:nrow(y)
    y = y[loc,]
    chrom = chrom[loc]
    
    T = nrow(y)

    
    if(!is.null(chrom)){
      temp=chrom[2:length(chrom)] - chrom[1:(length(chrom)-1)]
      chrom.start = which(temp==1)+1
      final = T
      chromnums = chrom[chrom.start]
      chrom.start = c(0, chrom.start, final)
      chromnums = c(chromnums[1]-1,chromnums)
    } else {
      chrom.start=vector(length=0)
    }
    
    if(horizontal){
      image.plot(loc,1:ncol(y),y,xlab=pos.string,ylab="Sample",...)

      if(length(chrom.start)>1){
        for(i in 1:(length(chrom.start)-1)){
          segments(chrom.start[i],0,chrom.start[i],nsamples)
          mtext(format(chromnums[i]),side=3,at=(chrom.start[i]+chrom.start[i+1])/2)
          if(i==length(chrom.start)-1){
            segments(chrom.start[length(chrom.start)],0,chrom.start[length(chrom.start)],nsamples)
          }
        }
      }
    } else {
      image.plot(1:ncol(y),loc,t(y),ylab=pos.string,xlab="Sample",...)
      if(length(chrom.start)>1){
        for(i in 1:(length(chrom.start)-1)){
          segments(0, chrom.start[i],nsamples, chrom.start[i])
          mtext(format(chromnums[i]),side=4,at=(chrom.start[i]+chrom.start[i+1])/2)
          if(i==length(chrom.start)-1){
            segments(0,chrom.start[length(chrom.start)],nsamples,chrom.start[length(chrom.start)])
          }
        }
      }
    }
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
