\name{compute.max.ProjectedZ}
\alias{compute.max.ProjectedZ}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ ~~function to do ... ~~ }
\description{
  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
compute.max.ProjectedZ(this.S, this.SST, this.imap, win, rratio, MIN.SNPs)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{this.S}{ ~~Describe \code{this.S} here~~ }
  \item{this.SST}{ ~~Describe \code{this.SST} here~~ }
  \item{this.imap}{ ~~Describe \code{this.imap} here~~ }
  \item{win}{ ~~Describe \code{win} here~~ }
  \item{rratio}{ ~~Describe \code{rratio} here~~ }
  \item{MIN.SNPs}{ ~~Describe \code{MIN.SNPs} here~~ }
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
function(this.S,this.SST,this.imap,win,rratio,MIN.SNPs){
    K=ncol(this.S) # number of platforms.
    
    if(is.na(rratio) || (length(rratio) != K)){
        rratio = rep(1,K)
    }
    
    
    this.Z = ComputeProjectedZ.fromS.R(this.S, this.SST, this.imap, win, rratio, MIN.SNPs)
    if(is.null(this.Z)){
       return(list(bestchpt=c(NA,NA), bestZ = NA, Z=NA))
        
    }
    
    this.Z = (this.Z - 1)/sqrt(2)    
 
    maxind = matrix.max(this.Z)
    bestchpt= c(maxind[1],maxind[1]+maxind[2])
    bestZ = this.Z[maxind[1],maxind[2]]
    
    list(bestchpt=bestchpt,bestZ=bestZ, Z=this.Z)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
