`mpcbs.mbic` <-
function(y,pos,anchor,MIN.SNPs=2, MIN.BP.LEN=1000, rratio=NULL, f=NULL,
                  MAX.CHPTS=30, platform.names=NULL,plots=TRUE,plotspdf=NULL,use.filtered.scan=TRUE){
                  
              #    MIN.SNPs=2; MIN.BP.LEN=1000; rratio=NULL; f=NULL; MAX.CHPTS=50; plots=TRUE; plotspdf=NULL; use.filtered.scan=TRUE
                  
    imap = anchor$imap
    
    cat("Using anchors containing ",length(anchor$merged.pos)," positions.\n",sep="")
    cat("Use filtered scan: ", use.filtered.scan,", maximum change-points: ", MAX.CHPTS,"\n")
    
    K=length(y)
    if(is.null(platform.names)){
        platform.names=vector("list",K)
        for(k in 1:K) platform.names[[k]] = paste("Platform",k)
    }
    
    if(is.null(rratio)) rratio = rep(1,K)

    T = length(anchor$merged.pos)
    win=T
    S = matrix(0,nrow=T,ncol=K)  # S[t,i] is the cumulative sum of platform i at position anchor$merged.pos[t]
    SST = rep(0, K)
    sigma = rep(0,K)
        
    for(i in 1:K){
        S.i = cumsum(y[[i]])
        S.i = c(0,S.i)
        S[,i] = S.i[imap[,i]] # S[t,i] sums up all probes in platform i falling before merged.pos[t] 
        SST[i] =  sum((y[[i]]-mean(y[[i]]))^2)
        sigma[i] = sqrt(compute.var(y[[i]]))
    }
    
    
    chpts = c(1,T)
    if(use.filtered.scan){
        bestZ=fcompute.max.ProjectedZ(S,SST,imap,win,rratio,MIN.SNPs,f=f)
    }else {
        bestZ =compute.max.ProjectedZ(S,SST,imap,win,rratio,MIN.SNPs)
    }
    
    best.subchpt= matrix(bestZ$bestchpt,ncol=1,nrow=2)
    best.Z = bestZ$bestZ
 
    splitnum=0
    chpt.hist = vector("list",MAX.CHPTS)
    mbic=vector("list",MAX.CHPTS)
    
    if(plots && !is.null(plotspdf)) pdf(plotspdf,height=10.5,width=8)
    
    while(length(chpts)<MAX.CHPTS+2){
           
        if(sum(!is.na(best.Z))==0) break
         
        # ---------------------------------------------------------------
        # Find the next split and compute left center, and right optimums.
        # ----------------------------------------------------------------
        splitnum=splitnum+1
        max.Z = max(best.Z,na.rm=TRUE)
        max.region = which.max(best.Z)
        newchpt = c(best.subchpt[1,max.region],best.subchpt[2,max.region])
        cat("Split ",splitnum,": ",newchpt[1],", ",newchpt[2],", Z-score = ",round(max.Z,1),".\n",sep="")
        
        S.L = S[chpts[max.region]:newchpt[1],]- matrix(S[chpts[max.region],],nrow=newchpt[1]-chpts[max.region]+1,ncol=K,byrow=T)
        S.M =  S[newchpt[1]:newchpt[2],]- matrix(S[newchpt[1],],nrow=newchpt[2]-newchpt[1]+1,ncol=K,byrow=T)
        S.R = S[newchpt[2]:chpts[max.region+1],] - matrix(S[newchpt[2],],nrow=chpts[max.region+1]-newchpt[2]+1,ncol=K,byrow=T)
        
        SST.L = rep(0,K); SST.M=rep(0,K); SST.R=rep(0,K)
        for(i in 1:K){
            SST.L[i] = sum( (y[[i]][imap[chpts[max.region],i]:(imap[newchpt[1],i]-1)] - mean(y[[i]][imap[chpts[max.region],i]:(imap[newchpt[1],i]-1)]))^2)
            SST.M[i] = sum( (y[[i]][imap[newchpt[1],i]:(imap[newchpt[2],i]-1)] - mean(y[[i]][imap[newchpt[1],i]:(imap[newchpt[2],i]-1)]))^2)
            SST.R[i] = sum( (y[[i]][imap[newchpt[2],i]:(imap[chpts[max.region+1],i]-1)] - mean(y[[i]][imap[newchpt[2],i]:(imap[chpts[max.region+1],i]-1)]))^2)
        }
        
        imap.L = imap[chpts[max.region]:newchpt[1],]
        imap.M = imap[newchpt[1]:newchpt[2],]
        imap.R = imap[newchpt[2]:chpts[max.region+1],]
        
        if(use.filtered.scan){
           bestZ.L = fcompute.max.ProjectedZ(S.L,SST.L,imap.L,win,rratio,MIN.SNPs,f=f)
           bestZ.M = fcompute.max.ProjectedZ(S.M,SST.M,imap.M,win,rratio,MIN.SNPs,f=f)
           bestZ.R = fcompute.max.ProjectedZ(S.R,SST.R,imap.R,win,rratio,MIN.SNPs,f=f)
        } else {
            bestZ.L = compute.max.ProjectedZ(S.L,SST.L,imap.L,win,rratio,MIN.SNPs)
            bestZ.M = compute.max.ProjectedZ(S.M,SST.M,imap.M,win,rratio,MIN.SNPs)
            bestZ.R = compute.max.ProjectedZ(S.R,SST.R,imap.R,win,rratio,MIN.SNPs)
        }
        
        # --------------------------------------
        # If the # base pairs in any of the new
        # segments is < MIN.BP.LEN set those Z
        # to NA.  don't need to change bestchpt
        # because those segments would never be selected.
        # this is different from when it's < MIN.SNPs,
        # in that case bestchpt would not be defined.
        # ---------------------------------------
        if(! is.na(bestZ.L$bestZ)){
            newlen.L = anchor$merged.pos[bestZ.L$bestchpt[1]+chpts[max.region]-1]-anchor$merged.pos[chpts[max.region]]
            newlen.M = anchor$merged.pos[bestZ.L$bestchpt[2]+chpts[max.region]-1]-anchor$merged.pos[bestZ.L$bestchpt[1]+chpts[max.region]-1] 
            newlen.R = anchor$merged.pos[newchpt[1]]-anchor$merged.pos[bestZ.L$bestchpt[2]+chpts[max.region]-1] 
            if( (newlen.M< MIN.BP.LEN)
            ||  ((newlen.L< MIN.BP.LEN) && (newlen.R< MIN.BP.LEN)) )
            {
                bestZ.L$bestZ = NA
            }
        }
        if(! is.na(bestZ.M$bestZ)){
            newlen.L = anchor$merged.pos[bestZ.M$bestchpt[1]+newchpt[1]-1]-anchor$merged.pos[newchpt[1]]
            newlen.M = anchor$merged.pos[bestZ.M$bestchpt[2]+newchpt[1]-1]-anchor$merged.pos[bestZ.M$bestchpt[1]+newchpt[1]-1] 
            newlen.R = anchor$merged.pos[newchpt[2]]-anchor$merged.pos[bestZ.M$bestchpt[2]+newchpt[1]-1] 
            if( (newlen.M< MIN.BP.LEN)
            ||  ((newlen.L< MIN.BP.LEN) && (newlen.R< MIN.BP.LEN)) )
            {
                bestZ.M$bestZ = NA
            }
        }
    
        if(! is.na(bestZ.R$bestZ)){
            newlen.L = anchor$merged.pos[bestZ.R$bestchpt[1]+newchpt[2]-1]-anchor$merged.pos[newchpt[2]]
            newlen.M = anchor$merged.pos[bestZ.R$bestchpt[2]+newchpt[2]-1]-anchor$merged.pos[bestZ.R$bestchpt[1]+newchpt[2]-1] 
            newlen.R = anchor$merged.pos[chpts[max.region+1]]-anchor$merged.pos[bestZ.R$bestchpt[2]+newchpt[2]-1] 
            if( (newlen.M< MIN.BP.LEN)
            ||  ((newlen.L< MIN.BP.LEN) && (newlen.R< MIN.BP.LEN) ) )
            {
                bestZ.R$bestZ = NA
            }
        }
        
        # -----------------------------------
        # Add new chpt into list.
        # -----------------------------------
        if(max.region>1){ 
            leftpart = best.subchpt[,1:(max.region-1)]
            leftpart.Z=best.Z[1:(max.region-1)]
        }else {
            leftpart = matrix(0,ncol=0, nrow=2)
            leftpart.Z = matrix(0,ncol=0,nrow=0)
        }
        if(max.region+1 <= ncol(best.subchpt)){ 
            rightpart = best.subchpt[,(max.region+1):ncol(best.subchpt)]
            rightpart.Z = best.Z[(max.region+1):length(best.Z)]
        }else  {
            rightpart = matrix(0,ncol=0, nrow=2)
            rightpart.Z=matrix(0,ncol=0,nrow=0)
        }
        chpt.hist[[splitnum]] = list(chpts=chpts,max.region=max.region, newchpt=newchpt,max.Z=max.Z)
        best.Z = c(leftpart.Z, bestZ.L$bestZ, bestZ.M$bestZ, bestZ.R$bestZ, rightpart.Z)    
        best.subchpt = cbind(leftpart, 
                             bestZ.L$bestchpt+chpts[max.region]-1, 
                             bestZ.M$bestchpt+newchpt[1]-1, 
                             bestZ.R$bestchpt+newchpt[2]-1,
                             rightpart)
        chpts = c(chpts[1:max.region],newchpt, chpts[(max.region+1):length(chpts)])  
        
        mbic[[splitnum]] = mbic.mp(S,imap=imap,sigma=sigma, tau=chpts,rratio=rratio)
    
        # ---------------------------------------
        # Progress plots.
        # ---------------------------------------
        if(plots){
      #      png(paste("mpcbs_",splitnum,".png",sep=""),height=1200,width=800)
            par(mfrow=c(K+1,1))
            notna = which(!is.na(best.Z))
            xmin=min(anchor$merged.pos); xmax=max(anchor$merged.pos)
            for(k in 1:K){
                plot(pos[[k]],y[[k]], pch=17,col="blue",xlim=c(xmin,xmax),main=platform.names[[k]],ylab="",xlab="base position")
                segments(anchor$merged.pos[chpts],rep(min(y[[k]]),length(chpts)),anchor$merged.pos[chpts],rep(max(y[[k]]),length(chpts)),col="black",lwd=2)
                if(length(notna)>0){
                    segments(anchor$merged.pos[best.subchpt[1,notna]],rep(min(y[[k]]),length(notna)),anchor$merged.pos[best.subchpt[1,notna]],rep(max(y[[k]]),length(notna)),col="red",lty=4)
                    segments(anchor$merged.pos[best.subchpt[2,notna]],rep(min(y[[k]]),length(notna)),anchor$merged.pos[best.subchpt[2,notna]],rep(max(y[[k]]),length(notna)),col="red",lty=4)
                }
            }
            Zplot = rep(0, length(anchor$merged.pos))
            if(length(notna)>0){
                for( i in 1:length(notna)){
                    Zplot[best.subchpt[1,notna[i]]:best.subchpt[2,notna[i]]] = sqrt(max(0,best.Z[notna[i]]))
                }
            }
            plot(anchor$merged.pos, Zplot, xlim=c(xmin,xmax),col="red",type="l",lwd=3,
                    main=paste("Cross-platform Evidence for the Next Split (Iteration ",splitnum,")",sep=""), ylab="- sqrt log p-value")
            segments(anchor$merged.pos[chpts],rep(min(Zplot),length(chpts)),anchor$merged.pos[chpts],rep(max(Zplot),length(chpts)),lwd=2,col="black")       
        }
        
        
    }
    
    # Summarize change-points.
    
    if(splitnum>0){
        mbic.tot = rep(0, splitnum)
        term1 = rep(0,splitnum)
        term2 = rep(0,splitnum)
        term3 = rep(0,splitnum)
        for(i in 1:splitnum){
            mbic.tot[i] = mbic[[i]]$mbic
            term1[i] = mbic[[i]]$term1
            term2[i] = mbic[[i]]$term2
            term3[i] = mbic[[i]]$term3
        }  
                
        kstar = which.max(mbic.tot)
        if(plots){
            par(mfrow=c(1,1))
            # plot(term1, type="b", ylim=c(min(c(term1,mbic.tot)), max(c(term1,mbic.tot))))
            plot(mbic.tot,type="b", col="red", main="Modified BIC")
            
            segments(kstar,0, kstar, max(mbic.tot))
            
            # plot(term2, type="b", col="blue",ylim=c(min(c(term2,term3)), max(c(term2,term3))))
            # points(term3,type="b", col="darkgreen")
        }

        
        chpts.final = c(chpt.hist[[kstar]]$chpts[1:chpt.hist[[kstar]]$max.region],
                        chpt.hist[[kstar]]$newchpt, 
                        chpt.hist[[kstar]]$chpts[(chpt.hist[[kstar]]$max.region+1):length(chpt.hist[[kstar]]$chpts)])  
    } else {
        chpt.final=c(1,T)
    }
    
    yhat = vector("list",K)
    for(k in 1:K){
        yhat[[k]]=rep(0,length(y[[k]]))
        yhat[[k]][1:imap[1,k]] =mean( y[[k]][1:(imap[chpts.final[2],k]-1)])
        for(i in 1:(length(chpts.final)-1)){
            yhat[[k]][imap[chpts.final[i],k]:(imap[chpts.final[i+1],k]-1)]=mean( y[[k]][imap[chpts.final[i],k]:(imap[chpts.final[i+1],k]-1)] )
        }
        
        yhat[[k]] = yhat[[k]][1:length(y[[k]])]  # sometimes, imap[chpts.final[i],k] = (total # snps in k platform)+1, even when i is not last change-point.
    }
    
    if(length(chpts.final)>2){
        chpts.final = chpts.final[2:(length(chpts.final)-1)]
    } else{
        chpts.final=integer(0)
    }
    
    if(plots) plot.crossplatform(pos,y,anchor,chpts.final, ranking=NULL, yhat=yhat,platform.names=platform.names)
    if(plots && !is.null(plotspdf)) dev.off()

    cat("MPCBS with MBIC regularization is done.\n")
    if(!is.null(plotspdf)) cat("  Progress plots are in the file ",plotspdf,".\n",sep="")
    
    # YS. Adding segment matrix to the output (with physical locations)
    physloc = anchor$merged.pos
    num.cuts = length(chpts.final)
    seg.bounds = chpts.final
    if(seg.bounds[1] != 1) seg.bounds = c(1,seg.bounds)
    if(seg.bounds[(num.cuts+1)] != length(physloc))  seg.bounds = c(seg.bounds,length(physloc))
    lenseg = length(seg.bounds)
    
    segment.mat = matrix(NA,nrow=(lenseg-1),ncol=2)        
    segment.mat[1,1] = 1
    segment.mat[,2] = seg.bounds[-1]
    fillCol = seg.bounds[-c(1,lenseg)]
    segment.mat[2:(lenseg-1),1] = fillCol+1        
       
    list(anchor=anchor, yhat=yhat,chpts=chpts.final, segmat=segment.mat, chpt.hist=chpt.hist, mbic=mbic.tot, term1=term1)
}

