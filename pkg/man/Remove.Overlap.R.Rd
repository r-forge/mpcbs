\name{Remove.Overlap.R}
\alias{Remove.Overlap.R}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ ~~function to do ... ~~ }
\description{
  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
Remove.Overlap.R(chpts)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
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
function(chpts){

    lookup <- nrow(chpts)
    
    ret.v <- rep(FALSE, lookup)
    ret.v[1] <- TRUE
    for(i in 2:lookup){
        overlap <- 0;
        st <- chpts[i,1];
        ed <- chpts[i,2];
        for(j in which(ret.v)){
            # check for overlap.
            if(     (st<=chpts[j,1] && ed>=chpts[j,2]) || 
                (st>=chpts[j,1] && st<=chpts[j,2]) || 
                (ed>=chpts[j,1] && ed<=chpts[j,2]) )
            {
                overlap <- 1;
                break;
            }
        }
        if(overlap == 0){
            ret.v[i] <- TRUE
        }
            
    }

    return(ret.v)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
