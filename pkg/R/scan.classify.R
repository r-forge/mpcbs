`scan.classify` <-
function(y,segs, CHISQ.PVAL.THRESH=0.001, MIN.REQ.ABSDIFF=0, MIN.SUFF.ABSDIFF=Inf){
    nsamples=dim(y)[2]
    nsnps=dim(y)[1]
    y.var <- compute.var(y)
    yhat = matrix(nrow=nrow(y), ncol=ncol(y),data=0)
    if(length(segs)>2){ throw = rep(FALSE,nrow(segs))
    }else  {throw = FALSE}
    
    
    if(length(segs)==2){
        y.class = scan.classify.seg(y,segs,y.var,CHISQ.PVAL.THRESH,MIN.REQ.ABSDIFF, MIN.SUFF.ABSDIFF)
        throw = sum(y.class)==0
    } else {
        y.class =matrix(ncol=nrow(segs),nrow=nsamples,data=FALSE)
        for(i in 1:nrow(segs)){
            seg=segs[i,]
            y.class[,i] <- scan.classify.seg(y,seg,y.var,CHISQ.PVAL.THRESH,MIN.REQ.ABSDIFF, MIN.SUFF.ABSDIFF)  #  This is the step that takes the most time.
            throw[i] = sum(y.class[,i])==0
            if(sum(y.class[,i])==1){
                if(seg[2]-seg[1]>1) {
                    yhat[(seg[1]+1):seg[2],y.class[,i]] = rep(mean(y[(seg[1]+1):seg[2],y.class[,i]]),seg[2]-seg[1])
                } else {
                    yhat[seg[1]+1,y.class[,i]] <- y[seg[1]+1,y.class[,i]]
                }
            } else {
                if(seg[2]-seg[1]>1) {
                    yhat[(seg[1]+1):seg[2],y.class[,i]] <- matrix(rep(apply(y[(seg[1]+1):seg[2],y.class[,i]],2,mean),seg[2]-seg[1]),nrow=seg[2]-seg[1],byrow=TRUE)
                } else {
                  yhat[seg[1]+1,y.class[,i]] <- y[seg[1]+1,y.class[,i]]
                }
            }
        }   
    }

    segs.f = segs[!throw,]
    y.class = y.class[,!throw]
    list(filtered.segs = segs.f, yhat = yhat, y.class = y.class)
}

