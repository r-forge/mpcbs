\name{scoreFit}
\alias{scoreFit}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ ~~function to do ... ~~ }
\description{
  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
scoreFit(segments, reps, trios, fraction.overlap = 0.8)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{segments}{ ~~Describe \code{segments} here~~ }
  \item{reps}{ ~~Describe \code{reps} here~~ }
  \item{trios}{ ~~Describe \code{trios} here~~ }
  \item{fraction.overlap}{ ~~Describe \code{fraction.overlap} here~~ }
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
function(segments,reps, trios, fraction.overlap=0.8){
    falsecalls=vector("list",length(segments));
    
    numsegparams=0;
    
    for( i in 1:length(segments)){
         segdim=ncol(segments[[i]])
         if( !is.null(segdim) ) numsegparams = segdim
         
    }
    if(numsegparams==0){
        cat("No calls were made.\n")
        score=0
    } else {
    
        # Replicates
    
        for (i in 1:nrow(reps)){
            seg1=segments[[reps[i,1]]]
            seg2=segments[[reps[i,2]]]
            
            if(length(seg1)==numsegparams) {
                nseg1 = 1
                seg1 = matrix(nrow=1, ncol=2, data=seg1)
            } else {
                nseg1=nrow(seg1)
            }
            if(length(seg2)==numsegparams) {
                nseg2 = 1
                seg2 = matrix(nrow=1, ncol=2, data=seg2)
            } else {
                nseg2=nrow(seg2)
            }
                
            if(nseg1 > 0 && nseg2>0){
                matchedpairs=matrix(nrow=nseg1,ncol=nseg2, data=0)
                for (j in 1:nseg1){
                    for (k in 1:nseg2){
                        overlap=intersect(c(seg1[j,1]:seg1[j,2]),c(seg2[k,1]:seg2[k,2]))
                        overlaplen=length(overlap)
                        if( (overlaplen/(seg1[j,2]-seg1[j,1]+1) >=fraction.overlap) && (overlaplen/(seg2[k,2]-seg2[k,1]+1) >=fraction.overlap)){
                            matchedpairs[j,k] = 1
                        }
                    }
                }
                # falsecalls{i} contains the segments for sample i that are not found
                # in its replicate. 
                # falsecalls{i} = ( segment parameters , P(FP+FN) )
                # In the case of replicates, if that segment is not found in the replicate
                # sample, then P(FP or FN) = 1, it is either a false positive in this
                # sample, or a false negative in the other sample.
            
                falsecalls[[reps[i,1]]]=seg1[which(apply(matchedpairs,1,sum)==0),]  
                if(length(falsecalls[[reps[i,1]]])==numsegparams) falsecalls[[reps[i,1]]] = matrix(nrow=1, ncol=2, data=falsecalls[[reps[i,1]]] )
                falsecalls[[reps[i,1]]]=cbind( falsecalls[[reps[i,1]]], rep(1,sum(apply(matchedpairs,1,sum)==0)))
                
                falsecalls[[reps[i,2]]]=seg2[which(apply(matchedpairs,2,sum)==0),] 
                if(length(falsecalls[[reps[i,2]]])==numsegparams) falsecalls[[reps[i,2]]] = matrix(nrow=1, ncol=2, data=falsecalls[[reps[i,2]]] )
            
                falsecalls[[reps[i,2]]]=cbind( falsecalls[[reps[i,2]]], rep(1,sum(apply(matchedpairs,2,sum)==0)))
            } 
            
            if(nseg1==0 && nseg2>0){
                falsecalls[[reps[i,2]]]=seg2
                falsecalls[[reps[i,2]]]=cbind( falsecalls[[reps[i,2]]], rep(1,nseg2))
            }
            if(nseg1>0 && nseg2==0){
                falsecalls[[reps[i,1]]]=seg1
                falsecalls[[reps[i,1]]]=cbind( falsecalls[[reps[i,1]]], rep(1,nseg1))
    
            }
            
        }
        
        # Trios
        
        for (i in 1:nrow(trios)){
            child=segments[[trios[i,1]]]
            p1=segments[[trios[i,2]]]
            p2=segments[[trios[i,3]]]
            
            if(length(child)==numsegparams) {
                nchild = 1
                child = matrix(nrow=1, ncol=2, data=child)
            } else {
                nchild=nrow(child)
            }
            if(length(p1)==numsegparams) {
                np1 = 1
                p1 = matrix(nrow=1, ncol=2, data=p1)
            } else {
                np1=nrow(p1)
            }
            if(length(p2)==numsegparams) {
                np2 = 1
                p2 = matrix(nrow=1, ncol=2, data=p2)
            } else {
                np2=nrow(p2)
            }

            
            # ------ calculate the child-p1, child-p2, p1-p2 matrices of matched aberations -----
            
            if(nchild>0){
                
                
                cp1matchedpairs=matrix(nrow=nchild,ncol=np1, data=0)
                if(np1>0){
                    for( j in 1:nchild ){
                        for( k in 1:np1 ){
                            overlap=intersect(c(child[j,1]:child[j,2]),c(p1[k,1]:p1[k,2]))
                            overlaplen=length(overlap)
                            if( (overlaplen/(child[j,2]-child[j,1]+1) >=fraction.overlap) && (overlaplen/(p1[k,2]-p1[k,1]+1) >=fraction.overlap) ){
                                cp1matchedpairs[j,k] = 1
                            }
                        }
                    }
                }
                
                
                cp2matchedpairs=matrix(nrow=nchild,ncol=np2, data=0)
                if(np2>0){
                    for( j in 1:nchild){
                        for( k in 1:np2){
                            overlap=intersect(c(child[j,1]:child[j,2]),c(p2[k,1]:p2[k,2]))
                            overlaplen=length(overlap)
                            if( (overlaplen/(child[j,2]-child[j,1]+1) >=fraction.overlap) && (overlaplen/(p2[k,2]-p2[k,1]+1) >=fraction.overlap)){
                                cp2matchedpairs[j,k] = 1
                            }
                        }
                    }
                }
                
                # ----- calculate the false calls in child, p1, and p2 -----
                # for each of c, p1, p2, compare with the other two members in the
                # family.
                
                childfc = rep(0,0)
                for (j in 1:nchild){
                    if ((np1==0 || sum(cp1matchedpairs[j,])==0) && (np2==0 || sum(cp2matchedpairs[j,])==0)){
                        # child T p1 F p2 F => P(FP+FN)=1
                        childfc=rbind(childfc, c(child[j,], 1))
                    }
                }
                falsecalls[[trios[i,1]]]=childfc
            
            }
            
        }
        
        score=0
        for( i in 1:length(falsecalls)){
            if (length(falsecalls[[i]])>0){
                temp=falsecalls[[i]]
                score=score+sum(temp[,numsegparams+1])
            }
        }

    }

    list(num.inconsistent=score, falsecalls=falsecalls)

  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
