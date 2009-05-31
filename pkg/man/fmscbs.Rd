\name{fmscbs}
\alias{fmscbs}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ ~~function to do ... ~~ }
\description{
  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
fmscbs(y, win, f = c(0.01, 0.1, 1), MIN.SNPs = 3, ALPHA = 0, GLOBAL.PVAL.CUTOFF = 1e-04, MAX.CHPTS = NA, WCHISQ.CUTOFF = NA, plots = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{y}{ ~~Describe \code{y} here~~ }
  \item{win}{ ~~Describe \code{win} here~~ }
  \item{f}{ ~~Describe \code{f} here~~ }
  \item{MIN.SNPs}{ ~~Describe \code{MIN.SNPs} here~~ }
  \item{ALPHA}{ ~~Describe \code{ALPHA} here~~ }
  \item{GLOBAL.PVAL.CUTOFF}{ ~~Describe \code{GLOBAL.PVAL.CUTOFF} here~~ }
  \item{MAX.CHPTS}{ ~~Describe \code{MAX.CHPTS} here~~ }
  \item{WCHISQ.CUTOFF}{ ~~Describe \code{WCHISQ.CUTOFF} here~~ }
  \item{plots}{ ~~Describe \code{plots} here~~ }
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
function(y,win,f=c(0.01,0.1,1), MIN.SNPs=3,ALPHA=0,GLOBAL.PVAL.CUTOFF=0.0001,MAX.CHPTS=NA,WCHISQ.CUTOFF=NA,plots=TRUE){
                                        # For debugging:
                                        # win=1000; Z=NULL; MIN.SNPs=3; ALPHA=0; MAX.CHPTS=NA; WCHISQ.CUTOFF=15; plots=TRUE; 

  N=dim(y)[2]
  T=dim(y)[1]
  DELTA = MIN.SNPs/T
  
  if(is.na(MAX.CHPTS)) MAX.CHPTS=floor(T/MIN.SNPs)
  
  if(is.na(WCHISQ.CUTOFF)){
    WCHISQ.CUTOFF = getCutoffMultisampleWeightedChisq(GLOBAL.PVAL.CUTOFF,T,DELTA,win,N,ALPHA)
    cat("MSCBS: weighted chisquare cutoff = ",WCHISQ.CUTOFF,"\n")
    
  }

  zlim=c(-1,1)

  y.var<-compute.var(y)
  yhat = matrix(rep(apply(y,2,mean),T),nrow=T,byrow=TRUE)
  y.r = y-yhat
  if(plots) image.plot(c(1:T),c(1:N),yhat,zlim=zlim,xlab="position",ylab="sample",main="Fitted");
  
  chpts = c(1,T)
#  bestZ =compute.max.Z.C(y.r,win,y.var,ALPHA,MIN.SNPs)
  bestZ =fcompute.max.Z(y.r,win,y.var,ALPHA,MIN.SNPs)
  best.subchpt= matrix(bestZ$bestchpt,ncol=1,nrow=2)
  best.Z = bestZ$bestZ


  splitnum=0
  chpt.hist = vector("list",MAX.CHPTS)
     
  while(TRUE){
    max.Z = max(best.Z,na.rm=TRUE)
    max.region = which.max(best.Z)
    
    if(max.Z<WCHISQ.CUTOFF) {
      cat("Maximum Z-score is ",max.Z,", which does not exceed cutoff of ",WCHISQ.CUTOFF,".  Segmentation finished.\n",sep="")
      break
    }

    if(length(chpts)>MAX.CHPTS+2){
      cat("Maximum number of change-points reached.  Segmentation finished.\n")
      break
    }
    if(is.na(best.subchpt[1,max.region])){
      cat("Optimal region has no valid change-points.  Segmentation finished.\n")      
      break
    }    
    splitnum=splitnum+1
    newchpt = c(best.subchpt[1,max.region],best.subchpt[2,max.region])

    # Classify samples and update yhat.
    y.r.classify = mscbs.classify(y.r[chpts[max.region]:chpts[(max.region+1)],] , newchpt-chpts[max.region]+1 , y.var , CHISQ.PVAL.THRESH=0.001)
    y.r.hat = matrix(0,nrow=T, ncol=N)
    y.r.hat[chpts[max.region]:chpts[(max.region+1)],] = y.r.classify$yhat
    yhat = yhat + y.r.hat
    y.r = y.r - y.r.hat
    if(plots) image.plot(c(1:T),c(1:N),yhat,zlim=zlim, xlab="position",ylab="sample",main="Fitted");
    
    
    if(!is.na(newchpt[2])){  # The added change consistes of two change-points.
      cat("Split ",splitnum,": ",newchpt[1],", ",newchpt[2],", Z-score = ",max.Z,".\n",sep="")
      y.r.L = y.r[chpts[max.region]:(newchpt[1]-1),]
      y.r.M =  y.r[newchpt[1]:(newchpt[2]-1),]
      y.r.R = y.r[newchpt[2]:chpts[max.region+1],]
      
#      bestZ.L = compute.max.Z.C(y.r.L,win,y.var,ALPHA,MIN.SNPs)
      bestZ.L = fcompute.max.Z(y.r.L,win,y.var,ALPHA,MIN.SNPs)
#      bestZ.M = compute.max.Z.C(y.r.M,win,y.var,ALPHA,MIN.SNPs)
      bestZ.M = fcompute.max.Z(y.r.M,win,y.var,ALPHA,MIN.SNPs)
#      bestZ.R = compute.max.Z.C(y.r.R,win,y.var,ALPHA,MIN.SNPs)
      bestZ.R = fcompute.max.Z(y.r.R,win,y.var,ALPHA,MIN.SNPs)
      best.Z.new=c(bestZ.L$bestZ, bestZ.M$bestZ, bestZ.R$bestZ)
      best.subchpt.new=cbind(bestZ.L$bestchpt+chpts[max.region]-1, 
        bestZ.M$bestchpt+newchpt[1]-1, 
        bestZ.R$bestchpt+newchpt[2]-1)
    } else { # The added change-point is a singleton.
      newchpt = newchpt[1]
      cat("Split ",splitnum,": ",newchpt,
          ", Z-score = ",max.Z,".\n",sep="")
      y.r.L = y.r[chpts[max.region]:(newchpt-1),]
      y.r.R = y.r[newchpt:chpts[max.region+1],]
      
      
#      bestZ.L = compute.max.Z.C(y.r.L,win,y.var,ALPHA,MIN.SNPs)
      bestZ.L = fcompute.max.Z(y.r.L,win,y.var,ALPHA,MIN.SNPs)
#      bestZ.R = compute.max.Z.C(y.r.R,win,y.var,ALPHA,MIN.SNPs)
      bestZ.R = fcompute.max.Z(y.r.R,win,y.var,ALPHA,MIN.SNPs)
      best.Z.new=c(bestZ.L$bestZ, bestZ.R$bestZ)
      best.subchpt.new=cbind(bestZ.L$bestchpt+chpts[max.region]-1, 
        bestZ.R$bestchpt+newchpt-1)
    }
    
    if(max.region>1){ 
      leftpart = best.subchpt[,1:(max.region-1)]
      leftpart.Z=best.Z[1:(max.region-1)]
    }else {
      leftpart = matrix(0,ncol=0, nrow=2)
      leftpart.Z = matrix(0,ncol=0,nrow=0)
    }
    if(max.region+1 <= ncol(best.subchpt)){ 
      rightpart = best.subchpt[,(max.region+1):ncol(best.subchpt)]
      rightpart.Z = best.Z[(max.region+1):length(best.Z)]
    }else  {
      rightpart = matrix(0,ncol=0, nrow=2)
      rightpart.Z=matrix(0,ncol=0,nrow=0)
    }
    chpt.hist[[splitnum]] = list(chpts=chpts,max.region=max.region,
               newchpt=newchpt,max.Z=max.Z,carriers=y.r.classify$carriers)
    best.Z = c(leftpart.Z, best.Z.new, rightpart.Z)    
    best.subchpt = cbind(leftpart,best.subchpt.new,rightpart)
    chpts = c(chpts[1:max.region],newchpt, chpts[(max.region+1):length(chpts)])

                                        # # For debugging: 
                                        #    cat(chpts)
                                        #    cat(rbind(best.Z, best.subchpt))
    
  }

  chpt.hist = chpt.hist[1:splitnum]
  if(length(chpts)>2){ 
    chpts = chpts[2:(length(chpts)-1)]
    
    # recompute yhat
    yhat[1:chpts[1],] = matrix(nrow=chpts[1], ncol=ncol(yhat),data=apply(as.matrix(y[1:chpts[1],]),2,mean), byrow=TRUE)
    if(length(chpts)>=2){
        for(i in 2:length(chpts)){
            yhat[(chpts[i-1]+1):chpts[i],] = matrix(nrow=chpts[i]-chpts[i-1], ncol=ncol(yhat),data=apply(as.matrix(y[(chpts[i-1]+1):chpts[i],]),2,mean), byrow=TRUE)    
        }
    }
    if(chpts[length(chpts)]<nrow(yhat)){
        yhat[(chpts[length(chpts)]+1):nrow(yhat),] = matrix(nrow=T-chpts[length(chpts)], ncol=ncol(yhat),data=apply(as.matrix(y[(chpts[length(chpts)]+1):nrow(yhat),]),2,mean), byrow=TRUE)    
    }
  } else{
     chpts = rep(0,0)
  }
  
  list(chpt.hist=chpt.hist, chpts=chpts,yhat=yhat)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
