`initialHyper.cnv2d` <-
function(cnvobj){
   if (!inherits(cnvobj, 'cnv2d')) 
      stop("First argument must be a cnv2d object.")
 
    p=0.00001
    a=0.999
    b=0.0005
    base = c(median(cnvobj$intensity[,1]+cnvobj$intensity[,2])/2,0)
    
    mu=c(0,0)
    v=diag(2)
    sigAA = diag(2)
    sigAB = diag(2)
    sigBA = diag(2)
    sigBB = diag(2)

    vvec = c(t(v))
    sigAAvec = c(t(sigAA))
    sigABvec = c(t(sigAB))
    sigBAvec = c(t(sigBA))    
    sigBBvec = c(t(sigBB))
    
    list(p=p,a=a,b=b,mu=mu,base=base,v=v,vvec=vvec,sigAA=sigAA,sigAB=sigAB,sigBA=sigBA,sigBB=sigBB,sigAAvec=sigAAvec,sigABvec=sigABvec,sigBAvec=sigBAvec,sigBBvec=sigBBvec)
}

