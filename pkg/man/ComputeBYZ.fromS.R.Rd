\name{ComputeBYZ.fromS.R}
\alias{ComputeBYZ.fromS.R}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ ~~function to do ... ~~ }
\description{
  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
ComputeBYZ.fromS.R(this.S, y.var, this.imap, win, delta, MIN.SNPs)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{this.S}{ ~~Describe \code{this.S} here~~ }
  \item{y.var}{ ~~Describe \code{y.var} here~~ }
  \item{this.imap}{ ~~Describe \code{this.imap} here~~ }
  \item{win}{ ~~Describe \code{win} here~~ }
  \item{delta}{ ~~Describe \code{delta} here~~ }
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
function(this.S,y.var,this.imap,win,delta,MIN.SNPs){
    T = nrow(this.S) # number of snps.
    N = ncol(this.S) # number of samples.
    dfnum <- 1;

    win = min(win, T-1)
    if(win==0){ 
        return(NULL)
    }

    l.delta.1=matrix(nrow=T, ncol=win, data=0);
    l.delta.2=matrix(nrow=T, ncol=win, data=0);
    Z=l.delta.1;
    for(i in 1:N){ 
        totalsnps = this.imap[T,i]-this.imap[1,i]
        dfden = totalsnps-2
        for(k in 1:win){
                nsnps = this.imap[(k+1):T,i]-this.imap[1:(T-k),i]
                diff1 <- this.S[(k+1):T, i]-this.S[1:(T-k),i]
                
                l.delta.1[1:(T-k),k]  = delta[1]*(diff1 - this.S[T,i]*nsnps/totalsnps)/sqrt(y.var[i]) - (delta[1]^2/2)*nsnps*(1-nsnps/totalsnps)
                l.delta.2[1:(T-k),k]  = delta[2]*(diff1 - this.S[T,i]*nsnps/totalsnps)/sqrt(y.var[i]) - (delta[2]^2/2)*nsnps*(1-nsnps/totalsnps)
                              
                set.to.zero = which(nsnps<MIN.SNPs | ((totalsnps-nsnps)<MIN.SNPs))
                l.delta.1[set.to.zero,k] = 0
                l.delta.2[set.to.zero,k] = 0
        }
        Z <- Z+pmax(l.delta.1,0)+pmax(l.delta.2,0)
    }
    return(Z)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
