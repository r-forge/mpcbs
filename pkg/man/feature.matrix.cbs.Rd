\name{feature.matrix.cbs}
\alias{feature.matrix.cbs}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ ~~function to do ... ~~ }
\description{
  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
feature.matrix.cbs(yhat, chpts)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{yhat}{ ~~Describe \code{yhat} here~~ }
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
function(yhat,chpts){
  T = nrow(yhat)
  N= ncol(yhat)
  if(chpts[1] != 1) chpts = c(1,chpts,T)

  feature.matrix = matrix(0,nrow=N,ncol=length(chpts)-1)
  regions = matrix(nrow=length(chpts)-1, ncol=2,data=0)
  for(i in 1:(length(chpts)-1)){
    feature.matrix[,i] = yhat[chpts[i],]
    regions[i,] = c(chpts[i], chpts[i+1])
  }
  
  list(feature.matrix=feature.matrix, feature.regions=regions)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
