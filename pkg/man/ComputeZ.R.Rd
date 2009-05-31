\name{ComputeZ.R}
\alias{ComputeZ.R}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ ~~function to do ... ~~ }
\description{
  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
ComputeZ.R(Y, T, win, ALPHA)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Y}{ ~~Describe \code{Y} here~~ }
  \item{T}{ ~~Describe \code{T} here~~ }
  \item{win}{ ~~Describe \code{win} here~~ }
  \item{ALPHA}{ ~~Describe \code{ALPHA} here~~ }
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
function(Y, T, win, ALPHA){
    dfnum <- 1;
    dfden <- T-2;



    log.alpha <- log(ALPHA)
    
    g <- function(u){
        i2 <- (u < 10 + log.alpha)
        r1 <- u
        r1[i2] <- r1[i2] * exp(u[i2]/2)/(ALPHA+exp(u[i2]/2))
        return( r1 )
    }
    

    N <- dim(y)[1]

    U=matrix(nrow=T, ncol=win, data=0);
    Z=U;


    for(i in 1:N){

        
        S=cumsum(y[i,]);
        SST = sum((y[i,]-S[T]/T)^2);
        
        for(k in 1:win){
            diff1 <- S[(k+1):T]-S[1:(T-k)]
                SSb = k*(diff1/k-S[T]/T)^2;
                SSb = SSb + (T-k)*((S[T]-diff1)/(T-k)-S[T]/T)^2;
                SSw = SST-SSb;
                U[1:(T-k),k] = (SSb/dfnum)/(SSw/dfden);
        }
        
        Z <- Z+g(U);  # Z[t,k]: change starting at t+1, ending at t+k.
    }
    return(Z)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
