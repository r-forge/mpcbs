\name{ComputeProjectedZ.fromS.R}
\alias{ComputeProjectedZ.fromS.R}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ ~~function to do ... ~~ }
\description{
  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
ComputeProjectedZ.fromS.R(this.S, this.SST, this.imap, win, rratio, MIN.SNPs)
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
    T = nrow(this.S) # Number of SNPs.
    N = ncol(this.S) # Number of samples/platforms

    win = min(win, T-1)
    if(win==0){ 
        return(NULL)
    }
    
    U=matrix(nrow=T, ncol=win, data=0);
    weightU = U;
    sum.weight.sq = 0
    Z=U;
    for(i in 1:N){
        totalsnps = this.imap[T,i]-this.imap[1,i]
        dfden = totalsnps-2
        for(k in 1:win){
                nsnps = this.imap[(k+1):T,i]-this.imap[1:(T-k),i]
                diff1 <- this.S[(k+1):T, i]-this.S[1:(T-k),i]
                
                temp<-diff1/nsnps-this.S[T,i]/totalsnps
                
                SSb = nsnps*(temp)^2;
                SSb = SSb + (totalsnps-nsnps)*((this.S[T,i]-diff1)/(totalsnps-nsnps)-this.S[T,i]/totalsnps)^2;
                SSw = this.SST[i]-SSb;
                U[1:(T-k),k] = sqrt(SSb/(SSw/dfden));
                weightU[1:(T-k),k] = sign(temp)*sqrt(nsnps*(1-nsnps/totalsnps)/(SSw/dfden))
#                weightU[1:(T-k),k] = sign(temp)
                
                set.to.zero = which(nsnps<MIN.SNPs | ((totalsnps-nsnps)<MIN.SNPs))
                U[set.to.zero,k] = 0
        }
        Z <- Z+rratio[i]*weightU*U;
        sum.weight.sq = sum.weight.sq + rratio[i]^2*weightU^2
    }
    Z<- Z^2/sum.weight.sq
    Z[which(is.na(Z), arr.ind=TRUE)]=0
    return(Z)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
