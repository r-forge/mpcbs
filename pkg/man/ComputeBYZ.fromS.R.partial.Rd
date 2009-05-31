\name{ComputeBYZ.fromS.R.partial}
\alias{ComputeBYZ.fromS.R.partial}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ ~~function to do ... ~~ }
\description{
  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
ComputeBYZ.fromS.R.partial(this.S, y.var, this.imap, start.inds, end.inds, delta, MIN.SNPs)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{this.S}{ ~~Describe \code{this.S} here~~ }
  \item{y.var}{ ~~Describe \code{y.var} here~~ }
  \item{this.imap}{ ~~Describe \code{this.imap} here~~ }
  \item{start.inds}{ ~~Describe \code{start.inds} here~~ }
  \item{end.inds}{ ~~Describe \code{end.inds} here~~ }
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
function(this.S,y.var,this.imap,start.inds, end.inds,delta,MIN.SNPs){
    T = nrow(this.S)
    N = ncol(this.S)
    dfnum <- 1;

    Z = matrix(nrow=length(start.inds),ncol=length(end.inds),data=0)

    totalsnps = this.imap[T,] - this.imap[1,]
    dfden = totalsnps-2
    for(i in 1:length(start.inds)){
      for(j in 1:length(end.inds)){
        st = start.inds[i]
        ed = end.inds[j]
        if(st > 0 && st < ed && ed<T){
          nsnps = this.imap[ed,] - this.imap[st,]
          diff1 = this.S[ed,]-this.S[st,]
          
          l.delta.1 = delta[1]*(diff1 - this.S[T,]*nsnps/totalsnps)/sqrt(y.var) - (delta[1]^2/2)*nsnps*(1-nsnps/totalsnps)
          l.delta.2 = delta[2]*(diff1 - this.S[T,]*nsnps/totalsnps)/sqrt(y.var) - (delta[2]^2/2)*nsnps*(1-nsnps/totalsnps)
          
          set.to.zero = which(nsnps<MIN.SNPs | ((totalsnps-nsnps)<MIN.SNPs))
          l.delta.1[set.to.zero] = 0
          l.delta.2[set.to.zero] = 0
          Z[i,j] = sum(pmax(l.delta.1,0)+pmax(l.delta.2,0))
        }
      }
    }

    Z
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
