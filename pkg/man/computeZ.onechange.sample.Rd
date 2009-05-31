\name{computeZ.onechange.sample}
\alias{computeZ.onechange.sample}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ ~~function to do ... ~~ }
\description{
  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
computeZ.onechange.sample(y, t, y.var)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{y}{ ~~Describe \code{y} here~~ }
  \item{t}{ ~~Describe \code{t} here~~ }
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
function(y,t,y.var){

  T =ncol(y)
#  cat("here!! T=", T,", t=", t,", nrow(y)=", nrow(y),"\n\n")

  St = apply(y[,1:t],1,sum)
  ST = apply(y,1,sum)
  Zsq = (St - (t/T)*ST)^2/(y.var*t*(1-t/T))
  pval = 1-pchisq(Zsq,1)

  list(Zsq, pval=pval)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
