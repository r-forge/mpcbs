`subsamp.errors` <-
function(errordata,nblocks,L){
    N=ncol(errordata$logR)
    M=nrow(errordata$logR)
    
    st = sample(M-L+1,nblocks,replace=TRUE)
    samps = sample(N,nblocks,replace=TRUE)
    
    temp.fun<-function(t){
        errordata$logR[t[1]:(t[1]+L-1),t[2]]
    }
    err.logR=apply(cbind(st,samps),1, "temp.fun")
    
    temp.fun<-function(t){
        errordata$bfreq[t[1]:(t[1]+L-1),t[2]]
    }
    err.bfreq=apply(cbind(st,samps),1, "temp.fun")
    
    temp.fun<-function(t){
        errordata$gtype[t[1]:(t[1]+L-1),t[2]]
    }
    err.gtype=apply(cbind(st,samps),1, "temp.fun")
     
    # String them together.
    err.logR.vec=matrix(err.logR,nrow=nblocks*L,ncol=1)
    err.bfreq.vec=matrix(err.bfreq,nrow=nblocks*L,ncol=1)
    err.gtype.vec=matrix(err.gtype,nrow=nblocks*L,ncol=1)
    
    list(logR=err.logR.vec, bfreq=err.bfreq.vec, gtype=err.gtype.vec)
}

