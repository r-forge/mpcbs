`compute.max.Z` <-
function(this.S,this.SST,this.imap,win,ALPHA,MIN.SNPs){
    K=ncol(this.S)
    
    this.Z = ComputeZ.fromS.R(this.S, this.SST, this.imap, win, ALPHA, MIN.SNPs)
    if(is.null(this.Z)){
       return(list(bestchpt=c(NA,NA), bestZ = NA, Z=NA))
        
    }
    
    g2<-function(u){ u^2*exp(u^2/2)/(ALPHA+exp(u^2/2))}
    g2.moments=computeMoments(g2,0);
    this.Z = (this.Z - K*g2.moments$psidot)/sqrt(g2.moments$psidotdot*K)    
 
    maxind = matrix.max(this.Z)
    bestchpt= c(maxind[1],maxind[1]+maxind[2])
    bestZ = this.Z[maxind[1],maxind[2]]
    
    list(bestchpt=bestchpt,bestZ=bestZ, Z=this.Z)
 
}

