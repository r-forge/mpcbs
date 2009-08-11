`readHyper` <-
function(filename){
    hypervec <- scan(filename)

    p = hypervec[1]
    a = hypervec[2]
    b = hypervec[3]
    mu = c(hypervec[4],hypervec[5])
    vvec = c(hypervec[6],hypervec[7],hypervec[8],hypervec[9])
    sigAAvec = c(hypervec[10],hypervec[11],hypervec[12],hypervec[13])
    sigABvec = c(hypervec[14],hypervec[15],hypervec[16],hypervec[17])
    sigBAvec = c(hypervec[18],hypervec[19],hypervec[20],hypervec[21])
    sigBBvec = c(hypervec[22],hypervec[23],hypervec[24],hypervec[25])
    base1 = hypervec[26]
    base2 = hypervec[27]
    
    v = matrix(vvec,nrow=2,ncol=2)
    sigAA = matrix(sigAAvec,nrow=2,ncol=2)
    sigAB = matrix(sigABvec,nrow=2,ncol=2)
    sigBA = matrix(sigBAvec,nrow=2,ncol=2)
    sigBB = matrix(sigBBvec,nrow=2,ncol=2)
    
    list(p=p,a=a,b=b,base=c(base1,base2), mu=mu,vvec=vvec,sigAAvec=sigAAvec,sigABvec=sigABvec,sigBAvec=sigBAvec,sigBBvec=sigBBvec,v=v,sigAA=sigAA,sigAB=sigAB,sigBA=sigBA,sigBB=sigBB)
}

