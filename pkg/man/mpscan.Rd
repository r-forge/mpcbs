\name{mpscan}
\alias{mpscan}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ ~~function to do ... ~~ }
\description{
  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
mpscan(y, pos, anchors, win, Z = NULL, MIN.SNPs = 1, ALPHA = 0, WCHISQ.CUTOFF = NA)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{y}{ ~~Describe \code{y} here~~ }
  \item{pos}{ ~~Describe \code{pos} here~~ }
  \item{anchors}{ ~~Describe \code{anchors} here~~ }
  \item{win}{ ~~Describe \code{win} here~~ }
  \item{Z}{ ~~Describe \code{Z} here~~ }
  \item{MIN.SNPs}{ ~~Describe \code{MIN.SNPs} here~~ }
  \item{ALPHA}{ ~~Describe \code{ALPHA} here~~ }
  \item{WCHISQ.CUTOFF}{ ~~Describe \code{WCHISQ.CUTOFF} here~~ }
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
function(y,pos,anchors,win,Z=NULL,MIN.SNPs=1, ALPHA=0,WCHISQ.CUTOFF=NA){
    anchor = merge.pos(pos,anchors=anchors)
    imap = anchor$imap
    
    K=length(y)

par(mfrow = c(K,1))
for(k in 1:K) plot(pos[[k]],y[[k]])


    T = length(anchor$merged.pos)
    S = matrix(0,nrow=T,ncol=K)  # S[t,i] is the cumulative sum of platform i at position anchor$merged.pos[t]
    SST = rep(0, K)
        
    
    for(i in 1:K){
        S.i = cumsum(y[[i]])
        S.i = c(0,S.i)
        S[,i] = S.i[imap[,i]]
        SST[i] =  sum((y[[i]]-mean(y[[i]]))^2);
plot(pos[[i]],S.i[2:length(S.i)],col="red")
lines(anchor$merged.pos,S[,i],col="green")   
    }
    
    # SST, S, index.map enable computation of the Z matrix.


    if(is.null(Z)){
        Z = ComputeZ.fromS.R(S, SST, imap, win, ALPHA, MIN.SNPs)
        g2<-function(u){ u^2*exp(u^2/2)/(ALPHA+exp(u^2/2))}
        g2.moments=computeMoments(g2,0);
        Z = (Z - K*g2.moments$psidot)/sqrt(g2.moments$psidotdot*K)    
    }

par(mfrow = c(K+1,1))
for(k in 1:K) plot(pos[[k]],y[[k]])
image.plot(anchor$merged.pos,1:win,sqrt(Z))


    Z.passinds = which(Z>WCHISQ.CUTOFF,arr.ind=TRUE)

    if(length(Z.passinds) == 0){
         y.hat <- y
         y.hat[,] <- mean(y)
         return( list(yhat=yhat, chpts=matrix(nrow=0, ncol=2), sampleswithsegment=NULL,chisq=NULL,Z=NULL)  )
    }

    Z.pass = Z[Z.passinds]
    Z.pass.order =order(Z.pass,decreasing=TRUE)
    Z.pass = Z.pass[Z.pass.order]
    Z.passinds = Z.passinds[Z.pass.order,]
    
    chpts = cbind(Z.passinds[,1]+1,Z.passinds[,1]+Z.passinds[,2])
    chpts.keep = Remove.Overlap.C(chpts)
    chpts = chpts[chpts.keep,]    

    if(length(chpts)==2) chpts=matrix(chpts,nrow=1)

par(mfrow = c(K,1))
for(k in 1:K){ 
    plot(pos[[k]],y[[k]])
    segments(anchor$merged.pos[chpts[,1]],rep(min(y[[k]]),nrow(chpts)),anchor$merged.pos[chpts[,1]],rep(max(y[[k]]),nrow(chpts)),col="blue")
    segments(anchor$merged.pos[chpts[,2]],rep(min(y[[k]]),nrow(chpts)),anchor$merged.pos[chpts[,2]],rep(max(y[[k]]),nrow(chpts)),col="red")
  }


    platform.chisq <- Chisq.Contrib.fromS.R(S, SST, imap, chpts) # chisq[i,j] holds contribution of platform j to segment i.

    yhat = vector("list",K)

    for(j in 1:K){
        yhat[[j]]=rep(0,length(y[[j]]))
        for(i in 1:nrow(chpts)){
            yhat[[j]][imap[chpts[i,1],j]:imap[chpts[i,2],j]] = mean(y[[j]][imap[chpts[i,1],j]:imap[chpts[i,2],j]])
        }
    }
    par(mfrow=c(K,1))
    for(j in 1:K){ 
        plot(pos[[j]],y[[j]])
        lines(pos[[j]],yhat[[j]],col="red",lwd=2)
    }
        
    list(yhat=yhat, chpts=chpts, platform.chisq = platform.chisq, Z=Z)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
