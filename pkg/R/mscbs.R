`mscbs` <-
function(y,win,MIN.SNPs=3,ALPHA=0,GLOBAL.PVAL.CUTOFF=0.0001,MAX.CHPTS=NA,WCHISQ.CUTOFF=NA,plots=TRUE){
                                        # For debugging:
                                        #  win=300; Z=NULL; MIN.SNPs=3; ALPHA=0; MAX.CHPTS=NA; WCHISQ.CUTOFF=NA; plots=TRUE; GLOBAL.PVAL.CUTOFF=0.0001

  N=dim(y)[2]
  T=dim(y)[1]
  DELTA = MIN.SNPs/T
  
  if(is.na(MAX.CHPTS)) MAX.CHPTS=floor(T/MIN.SNPs)
  
  if(is.na(WCHISQ.CUTOFF)){
    WCHISQ.CUTOFF = getCutoffMultisampleWeightedChisq(GLOBAL.PVAL.CUTOFF,T,DELTA,win,N,ALPHA)
    cat("MSCBS: weighted chisquare cutoff = ",WCHISQ.CUTOFF,"\n")
    
  }

  zlim=c(-3,3)

  y.var<-compute.var(y)
  yhat = matrix(rep(apply(y,2,mean),T),nrow=T,byrow=TRUE)
  y.r = y-yhat
  if(plots) image.plot(c(1:T),c(1:N),yhat,zlim=zlim,xlab="position",ylab="sample",main="Fitted");
  
  chpts = c(1,T)
  bestZ =compute.max.Z.C(y.r,win,y.var,ALPHA,MIN.SNPs)
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
    y.r.classify = mscbs.classify(y.r[chpts[max.region]:chpts[(max.region+1)],] , newchpt-chpts[max.region]+1 , y.var , CHISQ.PVAL.THRESH=0.001,)
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
      
      bestZ.L = compute.max.Z.C(y.r.L,win,y.var,ALPHA,MIN.SNPs)
      bestZ.M = compute.max.Z.C(y.r.M,win,y.var,ALPHA,MIN.SNPs)
      bestZ.R = compute.max.Z.C(y.r.R,win,y.var,ALPHA,MIN.SNPs)
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
      
      
      bestZ.L = compute.max.Z.C(y.r.L,win,y.var,ALPHA,MIN.SNPs)
      bestZ.R = compute.max.Z.C(y.r.R,win,y.var,ALPHA,MIN.SNPs)
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
  if(length(chpts)>2) chpts = chpts[2:length(chpts)-1]
  else chpts = rep(0,0)
  list(chpt.hist=chpt.hist, chpts=chpts,yhat=yhat)
}

