\name{scan.classify.seg}
\alias{scan.classify.seg}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ ~~function to do ... ~~ }
\description{
  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
scan.classify.seg(this.y, seg, y.var, CHISQ.PVAL.THRESH = 0.001, MIN.REQ.ABSDIFF = 0, MIN.SUFF.ABSDIFF = Inf)
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
function(this.y,seg,y.var,CHISQ.PVAL.THRESH=0.001, MIN.REQ.ABSDIFF=0, MIN.SUFF.ABSDIFF=Inf){
        
        sample.pval = computeZ.squarewave.sample(t(this.y),seg,y.var)$pval
        if(seg[2]-seg[1]>1){ 
            sample.med.diff=apply(as.matrix(this.y[(seg[1]+1):seg[2],]),2,median) - apply(as.matrix(this.y[-c((seg[1]+1):seg[2]),]),2,median)
        } else { 
            sample.med.diff=as.matrix(this.y[(seg[1]+1):seg[2],]) - apply(as.matrix(this.y[-c((seg[1]+1):seg[2]),]),2,median)
        }
       
        sample.abs.med.diff = abs(sample.med.diff)
        pass1 = sample.pval<CHISQ.PVAL.THRESH
        pass2 = sample.abs.med.diff > MIN.SUFF.ABSDIFF
        cut1 = sample.abs.med.diff  < MIN.REQ.ABSDIFF

        sampleswithsegment = (pass1 | pass2) & (!cut1)
        sampleswithsegment

  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
