\name{compute.var}
\alias{compute.var}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ ~~function to do ... ~~ }
\description{
  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
compute.var(y, use.mean = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{y}{ ~~Describe \code{y} here~~ }
  \item{use.mean}{ ~~Describe \code{use.mean} here~~ }
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
function(y, use.mean=TRUE){

    if(is.null(dim(y))){
         y = as.matrix(y,nrow=length(y),ncol=1)
    }

    N=dim(y)[2]
    T=dim(y)[1]
    y.diff<-y[1:(T-1),] - y[2:T,]
    
    if(!use.mean){
        if(N>1){
            y.var <- apply(y.diff^2,2,median)
        } else {
            y.var = median(y.diff^2)
        }
    } else {
        if(N>1){
            y.var <- apply(y.diff^2,2,mean)/2
        } else {
            y.var = mean(y.diff^2)/2
        }
    }
    
    y.var
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
