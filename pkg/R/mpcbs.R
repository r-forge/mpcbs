`mpcbs` <-
function(y,pos,anchors,win,Z=NULL,MIN.SNPs=2, MIN.BP.LEN=10000, ALPHA=0, rratio=NULL, f=NULL,
                  MAX.CHPTS=50, WCHISQ.CUTOFF=NA,platform.names=NULL,plots=TRUE,plotspdf=NULL,use.filtered.scan=TRUE, useProjectedChisquare=TRUE){
                                        # For debugging:
                                        #    anchors=c(1,3); win=500; Z=NULL; MIN.SNPs=2; ALPHA=0; MAX.CHPTS=20; WCHISQ.CUTOFF=20;
    
    anchor = merge.pos(pos,anchors=anchors)
    imap = anchor$imap
    
    cat("Created anchors containing ",length(anchor$merged.pos)," positions.\n",sep="")
    cat("Use filtered scan: ", use.filtered.scan,", use projected chisquare: ", useProjectedChisquare,"\n")
    
    K=length(y)
    if(is.null(platform.names)){
        platform.names=vector("list",K)
        for(k in 1:K) platform.names[[k]] = paste("Platform",k)
    }
    
    if(is.null(rratio)) rratio = rep(1,K)

    T = length(anchor$merged.pos)
    S = matrix(0,nrow=T,ncol=K)  # S[t,i] is the cumulative sum of platform i at position anchor$merged.pos[t]
    SST = rep(0, K)
        
    for(i in 1:K){
        S.i = cumsum(y[[i]])
        S.i = c(0,S.i)
        S[,i] = S.i[imap[,i]] 
        SST[i] =  sum((y[[i]]-mean(y[[i]]))^2);
    }
    
    
    chpts = c(1,T)
    if(useProjectedChisquare){
        if(use.filtered.scan){
            bestZ=fcompute.max.ProjectedZ(S,SST,imap,win,rratio,MIN.SNPs,f=f)
        }else {
            bestZ =compute.max.ProjectedZ(S,SST,imap,win,rratio,MIN.SNPs)
        }
    } else { 
        bestZ =compute.max.Z(S,SST,imap,win,ALPHA,MIN.SNPs)
    }
    best.subchpt= matrix(bestZ$bestchpt,ncol=1,nrow=2)
    best.Z = bestZ$bestZ
 
    splitnum=0
    chpt.hist = vector("list",MAX.CHPTS)
    
   if(plots && !is.null(plotspdf)) pdf(plotspdf,height=10.5,width=8)
    
    while(TRUE){
        max.Z = max(best.Z,na.rm=TRUE)
        max.region = which.max(best.Z)
        
        if(max.Z<WCHISQ.CUTOFF) break
        if(length(chpts) >MAX.CHPTS+2) break

        splitnum=splitnum+1
        
        newchpt = c(best.subchpt[1,max.region],best.subchpt[2,max.region])
        
        cat("Split ",splitnum,": ",newchpt[1],", ",newchpt[2],", Z-score = ",max.Z,".\n",sep="")
        
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
        
        if(useProjectedChisquare){
            if(use.filtered.scan){
               bestZ.L = fcompute.max.ProjectedZ(S.L,SST.L,imap.L,win,rratio,MIN.SNPs,f=f)
               bestZ.M = fcompute.max.ProjectedZ(S.M,SST.M,imap.M,win,rratio,MIN.SNPs,f=f)
               bestZ.R = fcompute.max.ProjectedZ(S.R,SST.R,imap.R,win,rratio,MIN.SNPs,f=f)
            } else {
                bestZ.L = compute.max.ProjectedZ(S.L,SST.L,imap.L,win,rratio,MIN.SNPs)
                bestZ.M = compute.max.ProjectedZ(S.M,SST.M,imap.M,win,rratio,MIN.SNPs)
                bestZ.R = compute.max.ProjectedZ(S.R,SST.R,imap.R,win,rratio,MIN.SNPs)
            }
        } else {
            bestZ.L = compute.max.Z(S.L,SST.L,imap.L,win,ALPHA,MIN.SNPs)
            bestZ.M = compute.max.Z(S.M,SST.M,imap.M,win,ALPHA,MIN.SNPs)
            bestZ.R = compute.max.Z(S.R,SST.R,imap.R,win,ALPHA,MIN.SNPs)
        }
        
        # Add new chpt into list.

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

        
        chpt.hist[[splitnum]] = list(chpts=chpts,max.region=max.region, newchpt=newchpt,max.Z=max.Z)
        best.Z = c(leftpart.Z, bestZ.L$bestZ, bestZ.M$bestZ, bestZ.R$bestZ, rightpart.Z)    
        best.subchpt = cbind(leftpart, 
                             bestZ.L$bestchpt+chpts[max.region]-1, 
                             bestZ.M$bestchpt+newchpt[1]-1, 
                             bestZ.R$bestchpt+newchpt[2]-1,
                             rightpart)
        chpts = c(chpts[1:max.region],newchpt, chpts[(max.region+1):length(chpts)])
    
    
        # plotting.
    
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
                    Zplot[best.subchpt[1,notna[i]]:best.subchpt[2,notna[i]]] = sqrt(best.Z[notna[i]])
                }
            }
            plot(anchor$merged.pos, Zplot, xlim=c(xmin,xmax),col="red",type="l",lwd=3,
                    main=paste("Cross-platform Evidence for the Next Split (Iteration ",splitnum,")",sep=""), ylab="- sqrt log p-value")
            segments(anchor$merged.pos[chpts],rep(min(Zplot),length(chpts)),anchor$merged.pos[chpts],rep(max(Zplot),length(chpts)),lwd=2,col="black")        
    
    #        par(mfrow=c(K,3))
    #        for(k in 1:K) {
    #             plot(S.L[,k], main=paste(bestZ.L$bestchpt[1],", ", bestZ.L$bestchpt[2],", ",bestZ.L$bestZ))
    #             segments(bestZ.L$bestchpt, rep(min(S.L[,k]),2), bestZ.L$bestchpt,rep(max(S.L[,k]),2),col="red")
    #             plot(S.M[,k], main=paste(bestZ.M$bestchpt[1],", ", bestZ.M$bestchpt[2],", ",bestZ.M$bestZ))
    #             segments(bestZ.M$bestchpt, rep(min(S.M[,k]),2), bestZ.M$bestchpt,rep(max(S.M[,k]),2),col="red")
    #             plot(S.R[,k], main=paste(bestZ.R$bestchpt[1],", ", bestZ.R$bestchpt[2],", ",bestZ.R$bestZ))
    #             segments(bestZ.R$bestchpt, rep(min(S.R[,k]),2), bestZ.R$bestchpt,rep(max(S.R[,k]),2),col="red")
    #        }
    
    #        dev.off()

        }
    }
    
    # Summarize change-points.
    
    chpts.final=c(1,T)
    ranking = c(0,0)
    maxZ = c(0,0)
    
    if(splitnum>0){
        for(i in 1:splitnum){
            mindist.1=min(abs(anchor$merged.pos[chpt.hist[[i]]$newchpt[1]]-anchor$merged.pos[chpts.final]))
            mindist.2=min(abs(anchor$merged.pos[chpt.hist[[i]]$newchpt[2]]-anchor$merged.pos[chpts.final]))
            if(mindist.1>MIN.BP.LEN) {
                chpts.final = c(chpts.final, chpt.hist[[i]]$newchpt[1])
                ranking = c(ranking, i)
                maxZ = c(maxZ, chpt.hist[[i]]$max.Z)
            }
            if(mindist.2>MIN.BP.LEN) {
                chpts.final = c(chpts.final, chpt.hist[[i]]$newchpt[2])
                ranking = c(ranking,i)
                maxZ = c(maxZ, chpt.hist[[i]]$max.Z)
            }
        }
    }
    
    ord = order(chpts.final)
    chpts.final = chpts.final[ord]
    ranking = ranking[ord]
    maxZ = maxZ[ord]
    
    yhat = vector("list",K)
    for(k in 1:K){
        yhat[[k]]=rep(0,length(y[[k]]))
        yhat[[k]][1:imap[1,k]] =mean( y[[k]][1:(imap[chpts.final[2],k]-1)])
        for(i in 1:(length(chpts.final)-1)){
            yhat[[k]][imap[chpts.final[i],k]:(imap[chpts.final[i+1],k]-1)]=mean( y[[k]][imap[chpts.final[i],k]:(imap[chpts.final[i+1],k]-1)] )
        }
    }
    
    if(length(chpts.final)>2){
        chpts.final = chpts.final[2:(length(chpts.final)-1)]
        ranking = ranking[2:(length(ranking)-1)]
        maxZ = maxZ[2:(length(maxZ)-1)]
    } else{
        chpts.final=integer(0)
        ranking=integer(0)
        maxZ=integer(0)
        
    }
    
    plot.crossplatform(pos,y,anchor,chpts.final, ranking,yhat,platform.names)
    
    if(plots && !is.null(plotspdf)) dev.off()

    cat("Segmentation finished.  Progress plots are in the file ",plotspdf,"\n",sep="")
    
    list(anchor=anchor,yhat=yhat,chpts=chpts.final,ranking=ranking,Z.score=maxZ)
}

