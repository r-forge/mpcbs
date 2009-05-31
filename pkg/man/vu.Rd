\name{vu}
\alias{vu}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ ~~function to do ... ~~ }
\description{
  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
vu(x, maxn = 1000, do.approx = (abs(x) < 0.2))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{ ~~Describe \code{x} here~~ }
  \item{maxn}{ ~~Describe \code{maxn} here~~ }
  \item{do.approx}{ ~~Describe \code{do.approx} here~~ }
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
function(x,maxn=1000,do.approx=(abs(x)<0.2)){
    x=as.matrix(x)
    vux = matrix(0,nrow=nrow(x),ncol=ncol(x))

#    if(sum(is.na(do.approx)) > 0 && msscan.debug.trace){
#        #check for errors here
#        print(do.approx)
#        scan()
#        print(x)
#        scan()
#    }


    if(is.logical(do.approx)){
        if(sum(do.approx) > 0)   vux[do.approx] = exp(-x[do.approx]*0.583);
    }else if(length(do.approx) > 0){
        vux[do.approx] = exp(-x[do.approx]*0.583);
    }
  

    if (sum(do.approx)<length(x)){
        notdo.approx.ix = which(!do.approx)
        n=matrix(c(1:maxn),nrow=1,ncol=maxn)
        summands = pnorm(-0.5*matrix(x[notdo.approx.ix],nrow=length(notdo.approx.ix),ncol=1)\%*\%sqrt(n))/matrix(rep(n,length(notdo.approx.ix)),ncol=length(n), byrow=TRUE)
        expterm = -2*apply(summands,1,"sum");
        vux[notdo.approx.ix] = (2/x[notdo.approx.ix]^2)*exp(expterm);
    }

    vux
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
