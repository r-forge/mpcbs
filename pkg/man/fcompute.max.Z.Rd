\name{fcompute.max.Z}
\alias{fcompute.max.Z}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ ~~function to do ... ~~ }
\description{
  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
fcompute.max.Z(this.y, win, y.var, ALPHA, MIN.SNPs, SINGLECHANGE.THRESH = 1e-04)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{this.y}{ ~~Describe \code{this.y} here~~ }
  \item{win}{ ~~Describe \code{win} here~~ }
  \item{y.var}{ ~~Describe \code{y.var} here~~ }
  \item{ALPHA}{ ~~Describe \code{ALPHA} here~~ }
  \item{MIN.SNPs}{ ~~Describe \code{MIN.SNPs} here~~ }
  \item{SINGLECHANGE.THRESH}{ ~~Describe \code{SINGLECHANGE.THRESH} here~~ }
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
function(this.y,win,y.var,ALPHA,MIN.SNPs,SINGLECHANGE.THRESH=0.0001){
    N=ncol(this.y)
    this.T=nrow(this.y)
    if(this.T<3*MIN.SNPs){
        
        this.Z = NA
        bestZ=NA
        bestchpt=c(NA,NA)

    } else {

        win = min(win, this.T-1)

        this.S = apply(this.y,2,cumsum)
        this.SST = apply(this.y,2,var)*(this.T-1)
        
        # Find [t1, t2] that maximizes Z using filtered scan.
        temp = fscan.max(this.S,this.SST,this.imap=NULL,MIN.SNPs=MIN.SNPs,ALPHA=ALPHA)
        bestchpt= temp$seg
        bestZ = temp$maxZ
        
        # Test the left change-point individually.
        if(bestZ!=0){
        
            if(bestchpt[1]==1){ 
                pval.L = 1
            } else {
                pval.L<-computeZ.onechange(t(this.y[1:bestchpt[2],]),bestchpt[1],y.var)$pval
            }
            if(bestchpt[2]==this.T){
                pval.R =1
            } else {
                pval.R<-computeZ.onechange(t(this.y[(bestchpt[1]+1):this.T,]),bestchpt[2]-bestchpt[1],y.var)$pval
            }
         # Pruning of left and right change-points (Olshen and Venkatraman suggestion) happens here, even though we haven't yet tested whether the double change-point has significant p-value.  The point is, in case the double change-point has significant p-value, then we would need to test for the significance of the left and right change-points anyway.
          if(pval.L>SINGLECHANGE.THRESH || bestchpt[1]<2*MIN.SNPs){
            # Added 6/15: we prune out change-points that are within MIN.SNPs to either end point.  In future, need
            # get rid of ##*## above and just use this step to control the boundary effects.
            
             #cat("Pruned out left change-point ", bestchpt[1]," out of ", this.T,"\n")
            bestchpt = c(bestchpt[2], NA)
          } else {
            if(pval.R>SINGLECHANGE.THRESH || (this.T-bestchpt[2]) < 2*MIN.SNPs){
#             cat("Pruned out right change-point ", bestchpt[2]," out of ", this.T,"\n")
             bestchpt = c(bestchpt[1], NA)
            }
          }
        } else {
          bestchpt = c(NA,NA)
        }
    }
    list(bestchpt=bestchpt,bestZ=bestZ)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
