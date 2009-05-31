\name{mscbs.classify}
\alias{mscbs.classify}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ ~~function to do ... ~~ }
\description{
  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
mscbs.classify(this.y, seg, y.var, CHISQ.PVAL.THRESH = 0.001, MIN.REQ.ABSDIFF = NA, MIN.SUFF.ABSDIFF = NA)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{this.y}{ ~~Describe \code{this.y} here~~ }
  \item{seg}{ ~~Describe \code{seg} here~~ }
  \item{y.var}{ ~~Describe \code{y.var} here~~ }
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
function(this.y,seg,y.var, CHISQ.PVAL.THRESH=0.001, MIN.REQ.ABSDIFF=NA, MIN.SUFF.ABSDIFF=NA){
#    cat("entered mscbs.classify!\n")
    
    N=dim(this.y)[2]
    T=dim(this.y)[1]

    if(is.na(MIN.REQ.ABSDIFF)){
      MIN.REQ.ABSDIFF = 0.5*sqrt(y.var)
    }
    if(is.na(MIN.SUFF.ABSDIFF)){
      MIN.SUFF.ABSDIFF = 5*sqrt(y.var)
    }
    
    if(is.na(seg[1]) && is.na(seg[2])){
      # Invalid change-point.
      sampleswithsegment = rep(0,N)
      this.yhat=matrix(0,nrow=nrow(this.y),ncol=ncol(this.y))
    } else {
      if(is.na(seg[2])){
        # Single change-point.
        # cat("Squeak!  ncol = ", ncol(this.y), " nrow=", nrow(this.y),"\n")
        sample.pval = computeZ.onechange.sample(t(this.y),seg[1],y.var)$pval
        pass1 = sample.pval<CHISQ.PVAL.THRESH
        pass2 = abs(apply(as.matrix(this.y[1:seg[1],]),2,mean) - apply(as.matrix(this.y[(seg[1]+1):T,]),2,mean))>MIN.SUFF.ABSDIFF
        cut1 = abs(apply(as.matrix(this.y[1:seg[1],]),2,mean) - apply(as.matrix(this.y[(seg[1]+1):T,]),2,mean))<MIN.REQ.ABSDIFF
        sampleswithsegment = (pass1 | pass2) & (!cut1)
        this.yhat=matrix(0,nrow=nrow(this.y),ncol=ncol(this.y))
        this.yhat[1:seg[1],sampleswithsegment] = matrix(rep(apply(as.matrix(this.y[1:seg[1],sampleswithsegment],nrow=seg[1]),2,mean),seg[1]),nrow=seg[1],byrow=TRUE)
        this.yhat[(seg[1]+1):T,sampleswithsegment] = matrix(rep(apply(as.matrix(this.y[(seg[1]+1):T,sampleswithsegment],nrow=T-seg[1]),2,mean),T-seg[1]),nrow=T-seg[1],byrow=TRUE)        

      } else {
        # Two change-points, square wave change, see who has change in middle.
        sample.pval = computeZ.squarewave.sample(t(this.y),seg,y.var)$pval
        pass1 = sample.pval<CHISQ.PVAL.THRESH
        pass2 = abs(apply(as.matrix(this.y[seg[1]:seg[2],]),2,mean) - apply(as.matrix(this.y[c(1:(seg[1]-1),(seg[2]+1):T),]),2,mean))>MIN.SUFF.ABSDIFF
        cut1 = abs(apply(as.matrix(this.y[seg[1]:seg[2],]),2,mean) - apply(as.matrix(this.y[c(1:(seg[1]-1),(seg[2]+1):T),]),2,mean))<MIN.REQ.ABSDIFF
        sampleswithsegment = (pass1 | pass2) & (!cut1)
        this.yhat=matrix(0,nrow=nrow(this.y),ncol=ncol(this.y))
        this.yhat[seg[1]:(seg[2]-1),sampleswithsegment] = matrix(rep(apply(as.matrix(this.y[seg[1]:(seg[2]-1),sampleswithsegment],nrow=seg[2]-seg[1]),2,mean),seg[2]-seg[1]),nrow=seg[2]-seg[1],byrow=TRUE)
        this.yhat[c(1:(seg[1]-1),(seg[2]):T),sampleswithsegment] = matrix(rep(apply(as.matrix(this.y[c(1:seg[1],(seg[2]):T),sampleswithsegment], nrow=T-(seg[2]-seg[1])),2,mean),T-(seg[2]-seg[1])),nrow=T-(seg[2]-seg[1]),byrow=TRUE)

      }
    }
    
#     cat("exiting mscbs.classify!\n")
   
    
    list(carriers=sampleswithsegment, yhat=this.yhat)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
