\name{norm.snp.pc}
\alias{norm.snp.pc}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ ~~function to do ... ~~ }
\description{
  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
norm.snp.pc(subdata, num.pc = 2, do.plots = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{subdata}{ ~~Describe \code{subdata} here~~ }
  \item{num.pc}{ ~~Describe \code{num.pc} here~~ }
  \item{do.plots}{ ~~Describe \code{do.plots} here~~ }
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
function(subdata, num.pc=2, do.plots=TRUE){

    nsamples = ncol(subdata)
    nsnps = nrow(subdata)
    subdata.norm = matrix(nrow=nrow(subdata),ncol=ncol(subdata),data=0)
    artifacts.all= matrix(nrow=nrow(subdata),ncol=ncol(subdata),data=0)
    
    # First, center subdata so that each SNP has median 0, and each sample has median 0.
    meds=apply(subdata,1,median)
    subdata = subdata - matrix(nrow=nrow(subdata),ncol=ncol(subdata),data=meds,byrow=FALSE)
    temp=apply(subdata,2,median)
    subdata = subdata - matrix(nrow=nrow(subdata),ncol=ncol(subdata),data=meds, byrow=TRUE)
    
    pc<-prcomp(t(subdata))

    if(do.plots){
        par(mfrow=c(num.pc,1))
        for(i in 1:num.pc){  plot(pc$rotation[,i]) }
    }

    for(k in 1:num.pc){
        artifacts.all = artifacts.all + t(pc$x[,k]\%*\%t(pc$rotation[,k]))
    }

    subdata.norm = subdata-artifacts.all

#    
#        T=1000
#        loc =(-T+1):0
#        
#        loc=loc+T
#        par(mfrow=c(3,1))
#        heatmap(subdata,loc=loc,zlim=c(-1,1))
#        heatmap(artifacts.all,loc=loc,zlim=c(-1,1))
#        heatmap(subdata.norm,loc=loc,zlim=c(-1,1))


    # Divide each SNP by its IQR.
    iqr = apply(subdata.norm, 1, IQR)
    iqr = iqr/median(iqr)
    iqr.mat = matrix(iqr, nrow=nrow(subdata),ncol=ncol(subdata), byrow=FALSE)
    subdata.norm = subdata.norm / iqr.mat
    
    
    list(normed=subdata.norm, artifacts = artifacts.all, snp.medians = meds, pc=pc)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
