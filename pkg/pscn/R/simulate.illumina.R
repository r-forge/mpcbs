`simulate.illumina` <-
function(T,mu.1, mu.2, errordata = NULL, L=100, R.BASELINE=2){
    
    nblocks = T/L
    colors=c("blue","red","darkgreen","black")
    softcolors=c("skyblue","pink","green","gray") 
    
    err = subsamp.errors(errordata,nblocks,L)
    
    # For those with undetermined genotype, just assign to closest cluster.
    nogtype=which(err$gtype==4)
    for(i in nogtype){
        d = err$bfreq[i]-c(0,0.5,1)
        if(sum(is.na(d))>0){
            err$gtype[i] = 1
        } else {
            err$gtype[i] = which.min(abs(d))
        }
    }
    
    # Transform to AB values.
    err.AB = illumina.getAB(err$logR,err$bfreq, R.BASELINE=R.BASELINE)

#    par(mfrow=c(2,2))
#    plot(err$logR)
#    plot(err$bfreq,col=colors[err$gtype])
#    plot(err.AB$A, col=colors[err$gtype], ylim=c(0,5))
#    plot(err.AB$B, col=colors[err$gtype], ylim=c(0,5))
    
    mean.A=(2-(err$gtype-1))*(R.BASELINE/2)
    mean.B=(err$gtype-1)*(R.BASELINE/2)
    res.A = err.AB$A - mean.A
    res.B = err.AB$B - mean.B
    
    # phased.gtype=1,2,3,4 for AA,AB,BA,BB
    phased.gtype = err$gtype
    phased.gtype[which(err$gtype==3)] = 4
    phased.gtype[which(err$gtype==2)] = sample(c(2,3),sum(err$gtype==2),replace=TRUE)
    
    mu.A = rep(0,T)
    mu.B = rep(0,T)
    grp1 = which(phased.gtype==1)
    mu.A[grp1] = mu.1[grp1]+mu.2[grp1]
    mu.B[grp1] = 0
    grp2 = which(phased.gtype==2)
    mu.A[grp2] = mu.1[grp2]
    mu.B[grp2] = mu.2[grp2]
    grp3 = which(phased.gtype==3)
    mu.A[grp3] = mu.2[grp3]
    mu.B[grp3] = mu.1[grp3]
    grp4 = which(phased.gtype==4)
    mu.A[grp4] = 0
    mu.B[grp4] = mu.1[grp4]+mu.2[grp4]
    
#    plot(mu.A,col=colors[phased.gtype])
#    plot(mu.B,col=colors[phased.gtype])
    
    A = res.A+mu.A
    B = res.B+mu.B
    
    
#    par(mfrow=c(2,2))
#    plot(err.AB$A, col=softcolors[phased.gtype], ylim=c(0,5),main="Subsampled A")
#    plot(err.AB$B, col=softcolors[phased.gtype], ylim=c(0,5),main="Subsampled B")
#    plot(A,col=softcolors[phased.gtype], ylim=c(0,5),main="A with added signal")
#    points(mu.A, pch="-", col=colors[phased.gtype])
#    plot(B,col=softcolors[phased.gtype], ylim=c(0,5),main="B with added signal")
#    points(mu.B, pch="-", col=colors[phased.gtype])

    temp=illumina.getlogRbfreq(A,B)
    
#    par(mfrow=c(2,2))
#    plot(err$logR)
#    plot(err$bfreq)
#    plot(temp$logR)
#    plot(temp$bfreq)
    
    list(A=A,B=B,phased.gtype = phased.gtype,logR=temp$logR,bfreq=temp$bfreq,mu.A=mu.A,mu.B=mu.B)
        
}

