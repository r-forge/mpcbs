\name{Chisq.Contrib.fromS.R}
\alias{Chisq.Contrib.fromS.R}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ ~~function to do ... ~~ }
\description{
  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
Chisq.Contrib.fromS.R(S, SST, imap, chpts)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{S}{ ~~Describe \code{S} here~~ }
  \item{SST}{ ~~Describe \code{SST} here~~ }
  \item{imap}{ ~~Describe \code{imap} here~~ }
  \item{chpts}{ ~~Describe \code{chpts} here~~ }
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
function(S,SST,imap,chpts){
    T = nrow(S)
    N = ncol(S)
    dfnum <- 1;
    n.segs <- nrow(chpts)
    
    ret.m <- matrix(nrow=n.segs, ncol=N, data=0)
 
    totalsnps = imap[T,]-imap[1,]
    
    for(i in 1:n.segs){
        st <- chpts[i,1]
        ed <- chpts[i,2]
        w=ed-st
        
        for(j in 1:N){
            nsnps = imap[ed,j]-imap[st,j]        
            mn1 <- S[T,j]/totalsnps[j]
            SSb <- (S[ed,j]-S[st,j] -nsnps*mn1)^2 / (nsnps*(1-nsnps/totalsnps[j]));
            SSw = SST[j]-SSb;
            ret.m[i,j] = (SSb/dfnum)/(SSw/(totalsnps[j]-2));
        }

    }
    
    return(ret.m)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
