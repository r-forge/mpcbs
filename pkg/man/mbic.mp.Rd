\name{mbic.mp}
\alias{mbic.mp}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ ~~function to do ... ~~ }
\description{
  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
mbic.mp(S, imap, sigma, tau, rratio)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{S}{ ~~Describe \code{S} here~~ }
  \item{imap}{ ~~Describe \code{imap} here~~ }
  \item{sigma}{ ~~Describe \code{sigma} here~~ }
  \item{tau}{ ~~Describe \code{tau} here~~ }
  \item{rratio}{ ~~Describe \code{rratio} here~~ }
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
function(S,imap,sigma,tau,rratio){

    m=length(tau)-2
    K=ncol(S)
    nrS = nrow(S)


    if(m>=1){
        # NOTE: this assumes that the last imap-position is at or to the right of the last SNP for all platforms.
        # this wouldn't be true if, e.g., only Agilent was used as probes, and the last position in anchor$merged.pos 
        # is still to the left of many of the Affy SNPs.  Then, the last row of imap would not contain the number of SNPs in Affy.
        tot.snps = imap[nrow(imap),]   # tot.snps here is actually tot.snps+1, verified that this is correct.
           
        X = matrix(ncol=K, nrow=m, data=0)
        n = matrix(ncol=K, nrow=m+1, data=0)
        delta = matrix(ncol=K, nrow=m, data=0)
        U = rep(0,m) 
        
        for(i in 1:m){
            for(k in 1:K){
                n[i,k] = imap[tau[i+1],k] 
            }
        }
        n[m+1,] = tot.snps
        
        for(i in 1:m){
            for(k in 1:K){
                X[i,k] = (S[tau[i+1],k] - n[i,k]*S[tau[i+2],k]/n[i+1,k])/(sigma[k]*sqrt(n[i,k]*(1-n[i,k]/n[i+1,k])))
                delta[i,k] = rratio[k]*sqrt(n[i,k])/sigma[k]
                if(!is.finite(X[i,k])){ 
                    X[i,k] = 0  # Some change-points may involve only points not in platform k, thus, the contribution from platform k should be 0.
                }
            }
            U[i] = sum(delta[i,]*X[i,])/sqrt(sum(delta[i,]^2))
        }
        
        effect.samp.size = rep(0,m+1)
        effect.samp.size[1] = sum(n[1,])
        for(i in 1:m){
            effect.samp.size[i+1] = sum(n[i+1,]-n[i,])
        }
        
        term1 = 0.5*sum(U^2)
        term2 = - 0.5*sum(log(effect.samp.size[effect.samp.size>0])) 
        term3 = - lchoose(nrS, m)
        
        mbic = term1+term2+term3
        
            # questions:  should there be a -log(tot.snps) term here?
            
    } else {
        mbic=0; term1 = 0;  term2 = 0; term3 = 0
    }
    
    list(mbic=mbic, term1=term1, term2=term2, term3=term3)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
