`vu` <-
function(x,maxn=1000,do.approx=(abs(x)<0.2)){
    x=as.matrix(x)
    vux = matrix(0,nrow=nrow(x),ncol=ncol(x))

#    if(sum(is.na(do.approx)) > 0 && msscan.debug.trace){
#        #check for errors here
#        print(do.approx)
#        scan()
#        print(x)
#        scan()
#    }


    if(is.logical(do.approx)){
        if(sum(do.approx) > 0)   vux[do.approx] = exp(-x[do.approx]*0.583);
    }else if(length(do.approx) > 0){
        vux[do.approx] = exp(-x[do.approx]*0.583);
    }
  

    if (sum(do.approx)<length(x)){
        notdo.approx.ix = which(!do.approx)
        n=matrix(c(1:maxn),nrow=1,ncol=maxn)
        summands = pnorm(-0.5*matrix(x[notdo.approx.ix],nrow=length(notdo.approx.ix),ncol=1)%*%sqrt(n))/matrix(rep(n,length(notdo.approx.ix)),ncol=length(n), byrow=TRUE)
        expterm = -2*apply(summands,1,"sum");
        vux[notdo.approx.ix] = (2/x[notdo.approx.ix]^2)*exp(expterm);
    }

    vux
}

