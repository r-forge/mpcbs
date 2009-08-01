`createHyper` <-
function(p, a, b, muvec, vvec, sigAAvec, sigABvec, sigBAvec, sigBBvec, basevec){
    v = matrix(vvec, nrow=2,ncol=2)
    sigAA = matrix(sigAAvec, nrow=2,ncol=2,byrow=TRUE)
    sigAB = matrix(sigABvec, nrow=2,ncol=2,byrow=TRUE)
    sigBA = matrix(sigBAvec, nrow=2,ncol=2,byrow=TRUE)
    sigBB = matrix(sigBBvec, nrow=2,ncol=2,byrow=TRUE)
    
    list(p=p,a=a,b=b,mu=muvec,base=basevec,v=v,vvec=vvec,sigAA=sigAA,sigAB=sigAB,sigBA=sigBA,sigBB=sigBB,sigAAvec=sigAAvec,sigABvec=sigABvec,sigBAvec=sigBAvec,sigBBvec=sigBBvec)
}

