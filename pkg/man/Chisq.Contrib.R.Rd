\name{Chisq.Contrib.R}
\alias{Chisq.Contrib.R}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ ~~function to do ... ~~ }
\description{
  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
Chisq.Contrib.R(y, T, chpts)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{y}{ ~~Describe \code{y} here~~ }
  \item{T}{ ~~Describe \code{T} here~~ }
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
function(y, T, chpts){
    dfnum <- 1;
    dfden <- T-2;

    N <- nrow(y)
    ret.m <- matrix(nrow=nrow(chpts), ncol=N, data=0)
    # chisq holds the contribution of each sample to the aberration.

    n.segs <- nrow(chpts)

    for(i in 1:n.segs){
        st <- chpts[i,1];
        ed <- chpts[i,2];
        w <- ed-st; # Changed 11/4, previously w<-ed-st+1.

        for(j in 1:N){
            mn1 <- sum(y[j,])/T
            SST <- sum( (y[j,]- mn1)^2 );
            SSb <- (( sum(y[j,(st+1):ed])-w*mn1)^2) / (w*(1-w/T));
            SSw = SST-SSb;
            ret.m[i,j] = (SSb/dfnum)/(SSw/dfden);
        }

    }
    return(ret.m)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
