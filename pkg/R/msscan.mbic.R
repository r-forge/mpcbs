`msscan.mbic` <-
function(y,y.var,win,ALPHA=0,verbose=TRUE,MIN.SNPS = 2,MAX.CHPTS=100,MAX.SKIPPED=MAX.CHPTS,do.plots=TRUE, PRIOR.CARRIER.PROB=0.1,PRIOR.SEGCOUNT=1,COMMON.PHAT=FALSE){
    
    y = t(y)
        
    N=dim(y)[1]
    T=dim(y)[2]
    if(verbose){
        cat("Running msscan.mbic: ",N," samples, ",T," SNPs, ALPHA=",ALPHA,"\n", sep="")
    }
    

    # ---- Compute Z matrix.  Z[t,k]: changed segment starting at t+1, ending at t+k. #
    Z = ComputeZ.C(y, T, win, ALPHA)
    if(MIN.SNPS>1) Z[,1:MIN.SNPS-1] = 0
    g2<-function(u){ u^2*exp(u^2/2)/(ALPHA+exp(u^2/2))}
    g2.moments=computeMoments(g2,0);
    Z = (Z - N*g2.moments$psidot)/sqrt(g2.moments$psidotdot*N)    
#    if(do.plots){    
#        par(mfrow=c(2,1))
#        heatmap(t(y))
#        image.plot(Z)
#    }
    Z.vec = matrix(data=t(Z),nrow=1)
    Z.vec.ind = matrix(data=c(1:prod(dim(Z))), nrow=nrow(Z),ncol=ncol(Z),byrow=TRUE)
            # check:  these two should be the same.
            #  Z.vec[Z.vec.ind[5,20:30]]
            #  Z[5,20:30]
    Z.vec.ord = order(Z.vec, decreasing=TRUE)
            # check: 
            # plot(Z.vec[Z.vec.ord])
    segbag = matrix(nrow=0,ncol=2)
    carriers = matrix(nrow=N,ncol=0)
    bic = rep(0,MAX.CHPTS)
    bic.all = matrix(nrow=N,ncol=MAX.CHPTS,data=0)
    sum.deltasq = rep(0,MAX.CHPTS)
    sum.deltasqwin = rep(0,MAX.CHPTS)
    Z.ind = 1
    Z.len = length(Z.vec)
    nskipped=0

    while(TRUE){
        # (1) Find the change-point with next largest Z that's not deleted in Z.vec.
        while(Z.ind<Z.len){
            if(Z.vec[Z.vec.ord[Z.ind]] != 0) break
            Z.ind = Z.ind+1
        }
        if(Z.ind >= Z.len || Z.vec[Z.vec.ord[Z.ind]]<0) break
        st = ceiling(Z.vec.ord[Z.ind]/win)
        w = Z.vec.ord[Z.ind]%%win
        if(w == 0) w = win
        ed = st + w
        st = st+1  # changed segment starts at st and ends at ed.
        if(verbose){
            cat("Considering segment: (",st,",",ed,")... ")
#            heatmap(t(y),loc=c((st-10):(ed+10)))
#            image.plot(Z[(st-10):(ed+10),])
        }
        
        chisq <- Chisq.Contrib.R(y, T, matrix(c(st-1,ed),nrow=1,ncol=2)) # chisq[i,j] holds contribution of sample j to segment i.
        samp.ord = order(chisq,decreasing=TRUE)  # the order that samples enter into bag.
        
        if(nrow(segbag)==0){
            bic.update = update.bic.baseline(t(y),segbag,carriers,rep(0,0),rep(0,0),c(st,ed),samp.ord,win,y.var, PRIOR.SEGCOUNT=5, do.plots=do.plots)
        } else {
            bic.update = update.bic.baseline(t(y),segbag,carriers,sum.deltasq[1:nrow(segbag)],
                                sum.deltasqwin[1:nrow(segbag)],c(st,ed),samp.ord,win,y.var,PRIOR.CARRIER.PROB=0.1,PRIOR.SEGCOUNT=1, do.plots=do.plots)
        }
        if(COMMON.PHAT){
            mbic.loc = bic.update$mbic.one.phat
        } else {    
            mbic.loc = bic.update$mbic.sep.phat
        }
        sum.deltasq.loc = bic.update$sum.deltasq.loc
        sum.deltasqwin.loc = bic.update$sum.deltasqwin.loc
        nhat = which.max(mbic.loc)
        
        if(mbic.loc[nhat]<0){
            # this change-point should not be added.
            # thus, overlaps with it should not be removed.
            if(verbose){
                cat("max bic over all subsets is negative.\n")
            }
            nskipped = nskipped+1;
        } else {
            # add change-point, update bic.
            bic.all[,nrow(segbag)+1] = mbic.loc
            bic[nrow(segbag)+1] = mbic.loc[nhat]
            sum.deltasq[nrow(segbag)+1] = sum.deltasq.loc[nhat]
            sum.deltasqwin[nrow(segbag)+1] = sum.deltasqwin.loc[nhat]
            
            # then, update segbag, carriers.
            segbag = rbind(segbag, c(st,ed))
            this.carriers= rep(0,N)
            this.carriers[samp.ord[1:nhat]] = 1
            carriers = cbind(carriers,this.carriers)
            
            # finally, remove overlaps from consideration.
            for(i in c(max((st-win+1),2):ed)){
                overlapping.w = max(1,st-i+1):win
                Z.vec[Z.vec.ind[i-1,overlapping.w]] = 0  # Because Z records segments in terms of  [st+1, win] 
            }
            if(verbose){
                cat("bic maxed at ",nhat," carriers, value = ",mbic.loc[nhat],".\n")
            }
        }
        if(verbose){
                cat("   ",nrow(segbag)," change-points in bag, ",nskipped," skipped.\n")
        }
        if(nrow(segbag)>=MAX.CHPTS  || nskipped>=MAX.SKIPPED) break
        Z.ind = Z.ind+1
    }
   
    # Choose the final number of change-points by maximizing MBIC, and populate yhat.
    
    yhat = matrix(nrow=nrow(y),ncol=ncol(y),data=apply(y,1,median)) 
    mhat = which.max(bic)
    if(bic[mhat]>0){
        # at least some CNVs are found, populate yhat.
        final.cnvs = segbag[1:mhat,]
        final.carriers = carriers[,1:mhat]
        for(i in 1:mhat){
            w=final.cnvs[i,2]-final.cnvs[i,1]+1
            if(w>1){
                if(sum(final.carriers[,i])>1){
                    yhat[which(final.carriers[,i]>0),final.cnvs[i,1]:final.cnvs[i,2]]= matrix(nrow=sum(final.carriers[,i]),ncol=w,data=apply(y[which(final.carriers[,i]>0),final.cnvs[i,1]:final.cnvs[i,2]],1,mean))
                } else {
                    yhat[which(final.carriers[,i]>0),final.cnvs[i,1]:final.cnvs[i,2]]= matrix(nrow=sum(final.carriers[,i]),ncol=w,data=mean(y[which(final.carriers[,i]>0),final.cnvs[i,1]:final.cnvs[i,2]]))
                }
            } else {
                yhat[which(final.carriers[,i]>0),final.cnvs[i,1]:final.cnvs[i,2]]= matrix(nrow=sum(final.carriers[,i]),ncol=w,data=y[which(final.carriers[,i]>0),final.cnvs[i,1]:final.cnvs[i,2]])          
            }
        }
    }
    
    if(do.plots){
        par(mfrow=c(2,1))
        plot(bic, type="b")
        lines(rep(mhat,2) ,c(min(bic),max(bic)), col="black")
        heatmap(t(carriers), main=paste("Carriers of CNVs in bag, first",mhat," selected."))
    }
    list(yhat=t(yhat), segbag=segbag, segbag.carriers=carriers, bic.all = bic.all, bic, bic=bic, chpts=final.cnvs, carriers=final.carriers)
}

