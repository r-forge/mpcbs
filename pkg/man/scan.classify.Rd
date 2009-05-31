\name{scan.classify}
\alias{scan.classify}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ ~~function to do ... ~~ }
\description{
  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
scan.classify(y, segs, CHISQ.PVAL.THRESH = 0.001, MIN.REQ.ABSDIFF = 0, MIN.SUFF.ABSDIFF = Inf)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{y}{ ~~Describe \code{y} here~~ }
  \item{segs}{ ~~Describe \code{segs} here~~ }
  \item{CHISQ.PVAL.THRESH}{ ~~Describe \code{CHISQ.PVAL.THRESH} here~~ }
  \item{MIN.REQ.ABSDIFF}{ ~~Describe \code{MIN.REQ.ABSDIFF} here~~ }
  \item{MIN.SUFF.ABSDIFF}{ ~~Describe \code{MIN.SUFF.ABSDIFF} here~~ }
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
function(y,segs, CHISQ.PVAL.THRESH=0.001, MIN.REQ.ABSDIFF=0, MIN.SUFF.ABSDIFF=Inf){
    nsamples=dim(y)[2]
    nsnps=dim(y)[1]
    y.var <- compute.var(y)
    yhat = matrix(nrow=nrow(y), ncol=ncol(y),data=0)
    if(length(segs)>2){ throw = rep(FALSE,nrow(segs))
    }else  {throw = FALSE}
    
    
    if(length(segs)==2){
        y.class = scan.classify.seg(y,segs,y.var,CHISQ.PVAL.THRESH,MIN.REQ.ABSDIFF, MIN.SUFF.ABSDIFF)
        throw = sum(y.class)==0
    } else {
        y.class =matrix(ncol=nrow(segs),nrow=nsamples,data=FALSE)
        for(i in 1:nrow(segs)){
            seg=segs[i,]
            y.class[,i] <- scan.classify.seg(y,seg,y.var,CHISQ.PVAL.THRESH,MIN.REQ.ABSDIFF, MIN.SUFF.ABSDIFF)  #  This is the step that takes the most time.
            throw[i] = sum(y.class[,i])==0
            if(sum(y.class[,i])==1){
                if(seg[2]-seg[1]>1) {
                    yhat[(seg[1]+1):seg[2],y.class[,i]] = rep(mean(y[(seg[1]+1):seg[2],y.class[,i]]),seg[2]-seg[1])
                } else {
                    yhat[seg[1]+1,y.class[,i]] <- y[seg[1]+1,y.class[,i]]
                }
            } else {
                if(seg[2]-seg[1]>1) {
                    yhat[(seg[1]+1):seg[2],y.class[,i]] <- matrix(rep(apply(y[(seg[1]+1):seg[2],y.class[,i]],2,mean),seg[2]-seg[1]),nrow=seg[2]-seg[1],byrow=TRUE)
                } else {
                  yhat[seg[1]+1,y.class[,i]] <- y[seg[1]+1,y.class[,i]]
                }
            }
        }   
    }

    segs.f = segs[!throw,]
    y.class = y.class[,!throw]
    list(filtered.segs = segs.f, yhat = yhat, y.class = y.class)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
