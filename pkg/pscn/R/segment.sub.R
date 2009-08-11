`segment.sub` <-
function(  cnvobj, verbose=TRUE, start.from.seg=FALSE, merge=TRUE, 
                    DALPHA.THRESH = 0.005, DBETA.THRESH = 0.005, DELTATHRESH = 5,   
                    BETA.THRESH=0.02,ALPHA.THRESH=0.02,MIN.SNPS=20,MIN.HET.SNPS=5,BASELINE.BETA.THRESH=0.1,FRAC.CHPTS=0.01){

# Hao 09/09/2008  DALPHA.THRESH and DBETA.THRESH changed from 0.01 to 0.005;  BETA.THRESH and ALPHA.THRESH changed from 0.05 to 0.02

# Hao 09/29/2008 MIN.SNPS changed from 20 to 50

# Hao 3/6/2009 MIN.SNPS changed from 50 to 20; MIN.HET.SNPS change from 10 to 5
      
      # merge=TRUE; start.from.seg=FALSE; DALPHA.THRESH = 0.005; DBETA.THRESH = 0.005; DELTATHRESH = 5;BETA.THRESH=0.02; ALPHA.THRESH=0.02; MIN.SNPS=10; MIN.HET.SNPS=5; BASELINE.BETA.THRESH=0.1; FRAC.CHPTS=0.01
                    
   if(is.null(cnvobj$sig)){
        cat("Unable to segment: need to run smooth.cnv2d first.\n")
        return
   }
   
   n=length(cnvobj)
   
   if(start.from.seg){
        alpha=cnvobj$theta.seg[,1]
        beta=cnvobj$theta.seg[,2]
   } else {
        alpha = cnvobj$theta[,1]
        beta = cnvobj$theta[,2]
   }
   dalpha = alpha[2:n] - alpha[1:n-1]
   dbeta = beta[2:n] - beta[1:n-1]
   estBS = cnvobj$estBS
   
   if(is.null(DALPHA.THRESH)) DALPHA.THRESH = quantile(abs(dalpha),1-FRAC.CHPTS)
   if(is.null(DBETA.THRESH)) DBETA.THRESH = quantile(abs(dbeta),1-FRAC.CHPTS)
   
#    xlim=c(0,length(cnvobj))
#    par(mfrow=c(2,1))
#    plot(alpha, main="alpha",xlim=xlim)
#    plot(dalpha, main="difference in alpha",xlim=xlim,ylim=c(-2*DALPHA.THRESH,2*DALPHA.THRESH))
#    abline(DALPHA.THRESH,0,col="red")
#    abline(-DALPHA.THRESH,0,col="red")
#    plot(beta,main="beta",xlim=xlim)
#    plot(dbeta,main="difference in beta",xlim=xlim,ylim=c(-2*DBETA.THRESH,2*DBETA.THRESH))
#    abline(DBETA.THRESH,0,col="red")
#    abline(-DBETA.THRESH,0,col="red")

    
   chpts = which(abs(dalpha)>DALPHA.THRESH | abs(dbeta) > DBETA.THRESH)+1
   nchpts = length(chpts)
   if(nchpts>1){
        segstart = rep(0,nchpts)  # where did this run of continuous changepoints start? (index in chpts)
        deltat = rep(0,nchpts)  
        segstart[1] = 1
        keepchpt = rep(0,nchpts)
        for(i in 2:nchpts){
                deltat[i] = chpts[i]-chpts[i-1]
                if (deltat[i]<=DELTATHRESH){
                    segstart[i] = segstart[i-1]
                } else {
                    segstart[i] = i
                }
        }
#       cbind(chpts,deltat,segstart)
        for(i in 2:nchpts){
            if( segstart[i]!= segstart[i-1]){
                # end of a run.
                keepchpt[i-1] = 1
                keepchpt[segstart[i-1]]=1
                keepchpt[i]=1
            }
        }
        finalchpts = chpts[which(keepchpt==1)]
   } else {
        finalchpts = chpts
   }
   
#    estCP = rep(0,n)
#    estCP[chpts]=1
#    estCP2 = rep(0,n)
#    estCP2[finalchpts] = 1
#    par(mfrow=c(4,1))
#    #plot(cnvobj,which.plot="allele1",loc=loc)
#    plot(cnvobj,which.plot="allele2",loc=loc,seg=FALSE)
##    plot(cnvobj,which.plot="alpha",loc=loc)
#    plot(cnvobj,which.plot="beta",loc=loc,seg=FALSE)
#    plot(estCP,xlim=c(min(loc),max(loc)))
#    plot(estCP2,xlim=c(min(loc),max(loc)))


    if(merge){
    
# Hao 09/09/2008 add 'if' and 'else'   
    if (length(finalchpts)>0){
        nfinal = length(finalchpts)
        segalpha = rep(0,nfinal+1)
        segalpha[1] = mean(cnvobj$theta[1:finalchpts[1],1])
        segbeta = rep(0,nfinal+1)
        segbeta[1] = mean(cnvobj$theta[1:finalchpts[1],2])
        seglen = rep(0,nfinal+1)
        seglen[1] = finalchpts[1]
        seghet = rep(0,nfinal+1)
        seghet[1] = sum(cnvobj$mixstate[1:finalchpts[1]]==2 | cnvobj$mixstate[1:finalchpts[1]]==3)
        if(nfinal>1){
            for(i in c(1:(nfinal-1))){
                segalpha[i+1] = mean(cnvobj$theta[finalchpts[i]:finalchpts[i+1],1])
                segbeta[i+1] = mean(cnvobj$theta[finalchpts[i]:finalchpts[i+1],2])
                seglen[i+1] = finalchpts[i+1]-finalchpts[i]+1
                seghet[i+1] = sum(cnvobj$mixstate[finalchpts[i]:finalchpts[i+1]]==2 | cnvobj$mixstate[finalchpts[i]:finalchpts[i+1]]==3)
      
            }
        }
        segalpha[nfinal+1] = mean(cnvobj$theta[finalchpts[nfinal]:length(cnvobj),1])
        segbeta[nfinal+1] = mean(cnvobj$theta[finalchpts[nfinal]:length(cnvobj),2])
        seglen[nfinal+1] = length(cnvobj) - finalchpts[nfinal]
        seghet[nfinal+1] = sum(cnvobj$mixstate[finalchpts[nfinal]:length(cnvobj)]==2 | cnvobj$mixstate[finalchpts[nfinal]:length(cnvobj)]==3)
      
        cbind(c(0,finalchpts),segalpha,segbeta,seglen,seghet)
        
      
        
        throwchpt = rep(0,nfinal)
        keepchpt = rep(0,nfinal) # can veto a throw.
        for(i in c(2:(nfinal+1))){
            dalpha = abs(segalpha[i]-segalpha[i-1])
            dbeta = abs(segbeta[i]-segbeta[i-1])
            
            if(dalpha<ALPHA.THRESH && dbeta<BETA.THRESH){
                # Not enough difference between the two segments.
                throwchpt[i-1] = 1
            } else { 
                if(seglen[i]<MIN.SNPS && i != nfinal+1){
                    dalpha2 = abs(segalpha[i+1]-segalpha[i-1])
                    dbeta2 = abs(segbeta[i+1]-segbeta[i-1])
                    if(dalpha2<ALPHA.THRESH && dbeta2<BETA.THRESH){
                        # Segment too small, and flanking regions have no change.  
                        # Throw out the entire segment.
                        throwchpt[i-1] = 1
                        throwchpt[i] = 1
                    } else{
                      # Flanking regions have changed, so this is a "fuzzy" transition due to low heterozygosity.
                        # Just throw out the change-point with lower probability.
                        if ( estBS[finalchpts[i-1]]<estBS[finalchpts[i]]) {
                            throwchpt[i]=1
                            keepchpt[i-1]=1
                        } else {
                            throwchpt[i-1]=1
                            keepchpt[i]=1
                        }
                    }
                } else if (dalpha<ALPHA.THRESH && seghet[i] <MIN.HET.SNPS){
                    # no change in ALPHA and there is not enough heterozygote snps to call a jump in beta.  
                    if(i<nfinal+1){
                        dalpha2 = abs(segalpha[i+1]-segalpha[i-1])
                        dbeta2 = abs(segbeta[i+1]-segbeta[i-1])
                        if(dalpha2<ALPHA.THRESH && dbeta2<BETA.THRESH){
                            # Flanking regions have no change.  
                            # Throw out the entire segment.
                            throwchpt[i-1] = 1
                            throwchpt[i] = 1
                        }  else {
                            # Flanking regions have changed, so this is a "fuzzy" transition due to low heterozygosity.
                            # Just throw out the change-point with lower probability.
                            if ( estBS[finalchpts[i-1]]<estBS[finalchpts[i]]) {
                                throwchpt[i]=1
                                keepchpt[i-1]=1
                            } else {
                                throwchpt[i-1]=1
                                keepchpt[i]=1
                            }
                        }
                    } else {
                        throwchpt[i-1]=1
                    }   
                }
            }
        }
 # Hao 09/16/2008 add two more criterions to throw points
       if (finalchpts[1]<=MIN.SNPS){
          throwchpt[1] = 1
          keepchpt[1] = 0
       }
       if (finalchpts[nfinal]>(length(cnvobj$intensity[,1])-MIN.SNPS+1)){
          throwchpt[nfinal] = 1
          keepchpt[nfinal] = 0
       }      
       
        
       if(verbose) cat("Throwing away ",sum(throwchpt & !keepchpt)," change-points due to outliers and local trends.\n")
        
        
  #      print(cbind(c(finalchpts,0),segalpha,segbeta,seglen,seghet,c(throwchpt,0)),digits=2)
        mergedchpts = finalchpts[which(throwchpt==0 | keepchpt)]
        
#     estCP2 = rep(0,n)
#    estCP2[mergedchpts] = 1
#    estCP3=rep(0,n)
#    estCP3[finalchpts]=1
#    #png(paste("sample",ind,"_",xlim[1],"_",xlim[2],".png",sep=""), height=1200, width=800)
#    par(mfrow=c(4,1))
#    plot(cnvobj,which.plot="allele1",loc=loc, seg=F)
#    #plot(cnvobj,which.plot="allele2",loc=loc)
# #   plot(cnvobj,which.plot="alpha",loc=loc)
#    plot(cnvobj,which.plot="beta",loc=loc,seg=F)
#    plot(estCP3,xlim=c(min(loc),max(loc)))
#    plot(estCP2,xlim=c(min(loc),max(loc)))
##    #dev.off()
        
        cnvobj$chpts = mergedchpts
    } else {
        mergedchpts = NULL 
      }       
        
    } else {
        cnvobj$chpts = finalchpts
    }

    
   # Hao 09/09/2008  add the 'if'  and  'else'
  if(length(mergedchpts)>0 && length(finalchpts)>0){ 
    
   chpts = c(1,cnvobj$chpts,length(cnvobj)+1)

   cnvobj$theta.seg = matrix(0,n,2)
   for(i in 2:length(chpts)){
    locs = (chpts[i-1]):(chpts[i]-1)
    cnvobj$theta.seg[locs,] = cbind(rep(mean(alpha[locs]),length(locs)), rep(mean(beta[locs]),length(locs)))
   }
   
   
  }else{
    cnvobj$chpts=NULL
    cnvobj$theta.seg = matrix(0,n,2)
    cnvobj$theta.seg[,2] = mean(beta)
    cnvobj$theta.seg[,1] = mean(alpha)
  }
  
   if(!is.null(BASELINE.BETA.THRESH)){
    baseline = which(abs(cnvobj$theta.seg[,2])<BASELINE.BETA.THRESH)
    cnvobj$theta.seg[baseline,2] = 0
   }
   
   
   cnvobj$sig.seg = matrix(0,n,2)
   mixstate=cnvobj$mixstate
   for(i in c(1:4)){
        cnvobj$sig.seg[mixstate==i,] = cnvobj$theta.seg[mixstate==i,]%*%cnvobj$statemat[[i]]
   }
   cnvobj
}

