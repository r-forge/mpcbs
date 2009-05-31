`norm.snp` <-
function(subdata, K=3){
    subdata=t(subdata)
    nsnps = ncol(subdata)
    
    
    # First, center subdata so that each SNP has median 0, and each sample has median 0.
    meds=apply(subdata,2,median)
    subdata = subdata - matrix(nrow=nrow(subdata),ncol=ncol(subdata),data=meds,byrow=TRUE)
    temp=apply(subdata,1,median)
    subdata = subdata - matrix(nrow=nrow(subdata),ncol=ncol(subdata),data=meds, byrow=FALSE)
    
    
    
    # Perform k-means clustering on the samples.
    clust = kmeans(subdata, centers=K)

#    clust = pam(subdata,k=K) # can also try kmedoids.
#    par(mfrow=c(1,1))
#    barplot(clust$cluster)
#    
#    T=2000
#    loc =(-T+1):0
#
#    loc=loc+T
#    heatmap(t(subdata[order(clust$cluster),loc]), zlim=c(-1,1))
#    for(k in 1:K) abline(sum(clust$cluster<=k)+0.5,0)
    
    # Within each cluster, normalize each SNP by its median.
    subdata.norm = matrix(0,nrow=nrow(subdata),ncol=nsnps)
    meds=matrix(0,nrow=K,ncol=nsnps)
    for(i in 1:K){
        inds = which(clust$cluster==i)
        # Subtract out median for each snp:
        meds[i,] = apply(subdata[inds,],2,median)
        meds.mat= matrix(meds[i,],nrow = length(inds),ncol=nsnps,byrow=TRUE)
        subdata.norm[inds,] = subdata[inds,] -meds.mat
        
#        par(mfrow=c(3,1))
#        loc=1:200
#        image(t(subdata[inds,loc]),zlim=c(-3,3))
#        image(t(meds.mat[,loc]),zlim=c(-3,3))
#        image(t(subdata.norm[inds,loc]),zlim=c(-3,3))
    }
    
    # Divide each SNP by its IQR.
    iqr = apply(subdata.norm, 2, IQR)
    iqr = iqr/median(iqr)
    iqr.mat = matrix(iqr, nrow=nrow(subdata),ncol=ncol(subdata), byrow=TRUE)
    subdata.norm = subdata.norm / iqr.mat
    
#    loc=1:200
#    par(mfrow=c(2,1))
#    image.plot(loc, 1:nrow(subdata),t(subdata[,loc]),xlab="SNP",ylab="Sample")
#    plot(loc,apply(subdata[,loc],2,median),type="l")
#    
#    par(mfrow=c(2,1))
#    colors=c("black","red","blue")
#    plot(meds[1,loc],type="l",col=colors[1])
#    for(i in 2:K) lines(meds[i,loc],col=colors[i])
#    plot(clust$centers[1,loc],type="l",col=colors[1])
#    for(i in 2:K) lines(clust$centers[i,loc],col=colors[i])
# Conclusion: not too much difference between median and mean.
    
#    snpid=125
#    hist(subdata[,snpid],20)



#    # look at the effects of normalization:
#       
#    T=2000
#    loc =(-T+1):0
#
#    loc=loc+T
#    par(mfrow=c(2,1))
#    heatmap(t(subdata[order(clust$cluster),loc]), zlim=c(-1,1))
#    for(k in 1:K) abline(sum(clust$cluster<=k)+0.5,0)
#    heatmap(t(subdata.norm[order(clust$cluster),loc]), zlim=c(-1,1))
#    for(k in 1:K) abline(sum(clust$cluster<=k)+0.5,0)


    list(normed=t(subdata.norm), cluster = clust$cluster, snp.medians = meds)
}

