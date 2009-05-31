`mscbs.classify` <-
function(this.y,seg,y.var, CHISQ.PVAL.THRESH=0.001, MIN.REQ.ABSDIFF=NA, MIN.SUFF.ABSDIFF=NA){
#    cat("entered mscbs.classify!\n")
    
    N=dim(this.y)[2]
    T=dim(this.y)[1]

    if(is.na(MIN.REQ.ABSDIFF)){
      MIN.REQ.ABSDIFF = 0.5*sqrt(y.var)
    }
    if(is.na(MIN.SUFF.ABSDIFF)){
      MIN.SUFF.ABSDIFF = 5*sqrt(y.var)
    }
    
    if(is.na(seg[1]) && is.na(seg[2])){
      # Invalid change-point.
      sampleswithsegment = rep(0,N)
      this.yhat=matrix(0,nrow=nrow(this.y),ncol=ncol(this.y))
    } else {
      if(is.na(seg[2])){
        # Single change-point.
        # cat("Squeak!  ncol = ", ncol(this.y), " nrow=", nrow(this.y),"\n")
        sample.pval = computeZ.onechange.sample(t(this.y),seg[1],y.var)$pval
        pass1 = sample.pval<CHISQ.PVAL.THRESH
        pass2 = abs(apply(as.matrix(this.y[1:seg[1],]),2,mean) - apply(as.matrix(this.y[(seg[1]+1):T,]),2,mean))>MIN.SUFF.ABSDIFF
        cut1 = abs(apply(as.matrix(this.y[1:seg[1],]),2,mean) - apply(as.matrix(this.y[(seg[1]+1):T,]),2,mean))<MIN.REQ.ABSDIFF
        sampleswithsegment = (pass1 | pass2) & (!cut1)
        this.yhat=matrix(0,nrow=nrow(this.y),ncol=ncol(this.y))
        this.yhat[1:seg[1],sampleswithsegment] = matrix(rep(apply(as.matrix(this.y[1:seg[1],sampleswithsegment],nrow=seg[1]),2,mean),seg[1]),nrow=seg[1],byrow=TRUE)
        this.yhat[(seg[1]+1):T,sampleswithsegment] = matrix(rep(apply(as.matrix(this.y[(seg[1]+1):T,sampleswithsegment],nrow=T-seg[1]),2,mean),T-seg[1]),nrow=T-seg[1],byrow=TRUE)        

      } else {
        # Two change-points, square wave change, see who has change in middle.
        sample.pval = computeZ.squarewave.sample(t(this.y),seg,y.var)$pval
        pass1 = sample.pval<CHISQ.PVAL.THRESH
        pass2 = abs(apply(as.matrix(this.y[seg[1]:seg[2],]),2,mean) - apply(as.matrix(this.y[c(1:(seg[1]-1),(seg[2]+1):T),]),2,mean))>MIN.SUFF.ABSDIFF
        cut1 = abs(apply(as.matrix(this.y[seg[1]:seg[2],]),2,mean) - apply(as.matrix(this.y[c(1:(seg[1]-1),(seg[2]+1):T),]),2,mean))<MIN.REQ.ABSDIFF
        sampleswithsegment = (pass1 | pass2) & (!cut1)
        this.yhat=matrix(0,nrow=nrow(this.y),ncol=ncol(this.y))
        this.yhat[seg[1]:(seg[2]-1),sampleswithsegment] = matrix(rep(apply(as.matrix(this.y[seg[1]:(seg[2]-1),sampleswithsegment],nrow=seg[2]-seg[1]),2,mean),seg[2]-seg[1]),nrow=seg[2]-seg[1],byrow=TRUE)
        this.yhat[c(1:(seg[1]-1),(seg[2]):T),sampleswithsegment] = matrix(rep(apply(as.matrix(this.y[c(1:seg[1],(seg[2]):T),sampleswithsegment], nrow=T-(seg[2]-seg[1])),2,mean),T-(seg[2]-seg[1])),nrow=T-(seg[2]-seg[1]),byrow=TRUE)

      }
    }
    
#     cat("exiting mscbs.classify!\n")
   
    
    list(carriers=sampleswithsegment, yhat=this.yhat)
}

