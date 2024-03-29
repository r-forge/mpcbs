`cross.platform.consensus` <-
function(segs, sigma2, platform.names, plots=TRUE){

    K=ncol(segs[[1]]$anchor$imap)
    
    yfit.all = matrix(nrow=0, ncol=K)
    colnames(yfit.all)<-platform.names
    seglen.all = matrix(nrow=0, ncol=K)
    colnames(seglen.all)<-platform.names

    seglabel = matrix(nrow=0, ncol=3)

    for(ind in c(1:length(segs))){
        chpts = c(1,segs[[ind]]$chpts)
        chpts.R = c(segs[[ind]]$chpts,  length(segs[[ind]]$anchor$merged.pos))
        yfit = matrix(nrow=length(chpts),ncol=K)
        seglen = matrix(nrow=length(chpts),ncol=K)
        
        for(i in 1:K){
            yhat = segs[[ind]]$yhat[[i]]
            yfit[,i] = yhat[segs[[ind]]$anchor$imap[chpts,i]]       
        }
        
        seglen = segs[[ind]]$anchor$imap[chpts.R,,drop=FALSE] - segs[[ind]]$anchor$imap[chpts,,drop=FALSE] 
        lab=cbind(rep(ind,length(length(chpts))), chpts, chpts.R)
        yfit.all = rbind(yfit.all, yfit)
        seglen.all = rbind(seglen.all, seglen)
        seglabel = rbind(seglabel, lab)
    }


    yfit = yfit.all
    seglen = seglen.all
    seglen[which(seglen==0, arr.ind=1)]= 0.01

    yfit[which(is.na(yfit),arr.ind=TRUE)]=0

    # compute the v matrix.
    v = matrix(data=sigma2, nrow=nrow(seglen), ncol=length(platform.names), byrow=TRUE)
    v = v/seglen
    sd = sqrt(v)

    # serial plot.
    
#    par(mfrow=c(1,1))
#    plot(yfit[,1],type="p",col="red", ylim=c(min(yfit, na.rm=TRUE),max(yfit, na.rm=TRUE)), xlab="Region index", ylab="Fitted values", pch=1, cex=1.5)
#    lines(yfit[,2],type="p",col="blue", pch=1, cex=1.5)
#    lines(yfit[,3],type="p",col="black", pch=1, cex=1.5)
#    for(i in 1:nrow(yfit)){
#        segments(i, yfit[i,1]-2*sd[i,1], i, yfit[i,1]+2*sd[i,1], col="red")
#        segments(i, yfit[i,2]-2*sd[i,2], i, yfit[i,2]+2*sd[i,2], col="blue")
#   #     segments(i, yfit[i,3]-2*sd[i,3], i, yfit[i,3]+2*sd[i,3], col="black")
#    }

    # Estimate the platform-specific signal ratio.
    if(K==1){
        cat("There is only one platform: r=1, there is nothing to estimate.\n")
    } else {
        cat("Estimating platform specific response ratio.\n")
        cat("For identifiability, the ratio for ",platform.names[[K]]," is fixed at 1.\n", sep="")
    }
    cpfit = wpca(yfit, v, x.lab=platform.names, plots=plots)
    
    consensus.cn = vector("list",length(segs))
    for(ind in 1:length(segs)){
        ccn = rep(0,length(segs[[ind]]$anchor$merged.pos))
        chpts = seglabel[seglabel[,1]==ind,]
        theta.fit = cpfit$u[seglabel[,1]==ind]
        if(nrow(chpts)>0){
            for(i in 1:nrow(chpts)){
                ccn[chpts[i,2]:chpts[i,3]] = theta.fit[i]+cpfit$alpha[K]
            }
        }
        consensus.cn[[ind]] =ccn
    }
    
    list(converged = cpfit$converged, consensus.cn=consensus.cn, alpha=cpfit$alpha, r=cpfit$r, fhat=yfit, v=v, theta=cpfit$u, fhat.fitted=cpfit$fitted)

}

