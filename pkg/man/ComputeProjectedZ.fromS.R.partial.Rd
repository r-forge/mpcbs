\name{ComputeProjectedZ.fromS.R.partial}
\alias{ComputeProjectedZ.fromS.R.partial}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ ~~function to do ... ~~ }
\description{
  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
ComputeProjectedZ.fromS.R.partial(this.S, this.SST, this.imap, start.inds, end.inds, rratio, MIN.SNPs)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{this.S}{ ~~Describe \code{this.S} here~~ }
  \item{this.SST}{ ~~Describe \code{this.SST} here~~ }
  \item{this.imap}{ ~~Describe \code{this.imap} here~~ }
  \item{start.inds}{ ~~Describe \code{start.inds} here~~ }
  \item{end.inds}{ ~~Describe \code{end.inds} here~~ }
  \item{rratio}{ ~~Describe \code{rratio} here~~ }
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
function(this.S,this.SST,this.imap,start.inds, end.inds,rratio,MIN.SNPs){
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
          temp = diff1/nsnps - this.S[T,]/totalsnps
          SSb = nsnps*(temp)^2
          SSb = SSb + (totalsnps-nsnps)*((this.S[T,]-diff1)/(totalsnps-nsnps)-this.S[T,]/totalsnps)^2;
          SSw = this.SST - SSb
          U = sqrt((SSb/dfnum)/(SSw/dfden))
          
          weightU = sign(temp)*sqrt(nsnps*(1-nsnps/totalsnps)/(SSw/dfden))
#          weightU = sign(temp)
          
          set.to.zero = which(nsnps<MIN.SNPs | ((totalsnps-nsnps)<MIN.SNPs))
          U[set.to.zero] = 0
          Z[i,j] = sum(U*weightU*rratio)/sqrt(sum(rratio^2*weightU^2))
        }
      }
    }
    
    Z<- Z^2
    Z[which(is.na(Z), arr.ind=TRUE)]=0
 
    Z
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
