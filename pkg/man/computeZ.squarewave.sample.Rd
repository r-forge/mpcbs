\name{computeZ.squarewave.sample}
\alias{computeZ.squarewave.sample}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ ~~function to do ... ~~ }
\description{
  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
computeZ.squarewave.sample(this.y, seg, y.var)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{this.y}{ ~~Describe \code{this.y} here~~ }
  \item{seg}{ ~~Describe \code{seg} here~~ }
  \item{y.var}{ ~~Describe \code{y.var} here~~ }
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
function(this.y,seg,y.var){
  T =ncol(this.y)
  Ss = apply(as.matrix(this.y[,1:seg[1]],nrow=nrow(this.y),ncol=seg[1]),1,sum)
  St = apply(this.y[,1:seg[2]],1,sum)
  ST = apply(this.y,1,sum)
  k=seg[2]-seg[1]
  Zsq = (St-Ss - (k/T)*ST)^2/(y.var*k*(1-k/T))
  pval = 1-pchisq(Zsq,1)

  list(Zsq=Zsq, pval=pval)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
