\name{cnvtable}
\alias{cnvtable}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ ~~function to do ... ~~ }
\description{
  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
cnvtable(chrom, pos, featuremat, segments, sampleswithsegment)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{chrom}{ ~~Describe \code{chrom} here~~ }
  \item{pos}{ ~~Describe \code{pos} here~~ }
  \item{featuremat}{ ~~Describe \code{featuremat} here~~ }
  \item{segments}{ ~~Describe \code{segments} here~~ }
  \item{sampleswithsegment}{ ~~Describe \code{sampleswithsegment} here~~ }
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
function(chrom,pos,featuremat,segments,sampleswithsegment){
    nsegs = nrow(segments)
    
    NCOLS = 11
    cnvlist=matrix(0,nrow=0,ncol=NCOLS)
   
    for(i in 1:nsegs){
        st = segments[i,1]; ed = segments[i,2]
        carriers = which(sampleswithsegment[,i]==1)
        n.with.cnv = length(carriers)
        carriers.gain = carriers[which(featuremat[carriers,i]>0)]
        carriers.loss = carriers[which(featuremat[carriers,i]<0)]
        if(length(carriers.gain)>0) mean.gain = mean(featuremat[carriers.gain,i])
        else mean.gain=0
        if(length(carriers.loss)>0) mean.loss = mean(featuremat[carriers.loss,i])
        else mean.loss=0
        cnvlist = rbind(cnvlist,c(st,ed,ed-st+1,chrom[st],pos[st],pos[ed],length(carriers),length(carriers.gain),mean.gain,length(carriers.loss),mean.loss))
    }
    cnvlist=as.data.frame(cnvlist)
    names(cnvlist)<-c("Start index","End index","# SNPs","Chromosome","Start pos","End pos","# Carriers","# Gain","Mean gain","# Loss","Mean loss")
  
    cnvlist
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
