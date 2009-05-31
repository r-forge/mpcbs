`norm.snp.pc` <-
function(subdata, num.pc=2, do.plots=TRUE){

    nsamples = ncol(subdata)
    nsnps = nrow(subdata)
    subdata.norm = matrix(nrow=nrow(subdata),ncol=ncol(subdata),data=0)
    artifacts.all= matrix(nrow=nrow(subdata),ncol=ncol(subdata),data=0)
    
    # First, center subdata so that each SNP has median 0, and each sample has median 0.
    meds=apply(subdata,1,median)
    subdata = subdata - matrix(nrow=nrow(subdata),ncol=ncol(subdata),data=meds,byrow=FALSE)
    temp=apply(subdata,2,median)
    subdata = subdata - matrix(nrow=nrow(subdata),ncol=ncol(subdata),data=meds, byrow=TRUE)
    
    pc<-prcomp(t(subdata))

    if(do.plots){
        par(mfrow=c(num.pc,1))
        for(i in 1:num.pc){  plot(pc$rotation[,i]) }
    }

    for(k in 1:num.pc){
        artifacts.all = artifacts.all + t(pc$x[,k]%*%t(pc$rotation[,k]))
    }

    subdata.norm = subdata-artifacts.all

#    
#        T=1000
#        loc =(-T+1):0
#        
#        loc=loc+T
#        par(mfrow=c(3,1))
#        heatmap(subdata,loc=loc,zlim=c(-1,1))
#        heatmap(artifacts.all,loc=loc,zlim=c(-1,1))
#        heatmap(subdata.norm,loc=loc,zlim=c(-1,1))


    # Divide each SNP by its IQR.
    iqr = apply(subdata.norm, 1, IQR)
    iqr = iqr/median(iqr)
    iqr.mat = matrix(iqr, nrow=nrow(subdata),ncol=ncol(subdata), byrow=FALSE)
    subdata.norm = subdata.norm / iqr.mat
    
    
    list(normed=subdata.norm, artifacts = artifacts.all, snp.medians = meds, pc=pc)
}

