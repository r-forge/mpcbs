\name{merge.pos}
\alias{merge.pos}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ ~~function to do ... ~~ }
\description{
  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
merge.pos(pos, anchors = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{pos}{ ~~Describe \code{pos} here~~ }
  \item{anchors}{ ~~Describe \code{anchors} here~~ }
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
function(pos,anchors=NULL){
    M = length(pos)
    pos.len = rep(0,M)
    for(i in 1:M) pos.len[i] = length(pos[[i]])
    
    if (is.null(anchors)) anchors=1:M
    
    nanchors = length(anchors)
    not.anchors = setdiff(c(1:M),anchors)
    
    totalinds = 0
    for(i in 1:nanchors){
        totalinds =totalinds + length(pos[[anchors[i]]]) 
    }
    
    merged.pos = rep(0,totalinds)
    index.map = matrix(0,nrow=totalinds,ncol=M)
    
    counter = rep(1,M) # counter[i] keeps track of where current merged.inds is in each of pos[i].
    next.pos = rep(0,nanchors)
    for(i in 1:nanchors) next.pos[i] = pos[[anchors[i]]][counter[anchors[i]]]
    
    for(k in 1:totalinds){
        min.pos = which.min(next.pos)
        merged.pos[k] = next.pos[min.pos]
        next.pos[min.pos] = pos[[anchors[min.pos]]][counter[anchors[min.pos]]+1]
        counter[anchors[min.pos]] = counter[anchors[min.pos]]+1
        for(j in not.anchors){
            while(TRUE) {
                if(counter[j] > pos.len[j]) break
                if(pos[[j]][counter[j]]>merged.pos[k]) break
                counter[j] = counter[j] + 1
            }
        }
        
        index.map[k,] = counter
    }

    # merged.pos is union of pos[[anchors]].
    # imap[k,i] = min_j { pos[[i]][j] > merged.pos[k] }   (note it is STRICTLY larger than)  
    # Thus to get the index of the largest value in pos[[i]] 
    # that is smaller than or equal to merged.pos[k], use imap[i,k]-1.
    # By this notation, imap[totalinds,] = pos.len + 1.

    list(merged.pos = merged.pos, imap = index.map)

  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
