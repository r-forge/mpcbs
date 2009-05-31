`compute.max.ProjectedZ` <-
function(this.S,this.SST,this.imap,win,rratio,MIN.SNPs){
    K=ncol(this.S) # number of platforms.
    
    if(is.na(rratio) || (length(rratio) != K)){
        rratio = rep(1,K)
    }
    
    
    this.Z = ComputeProjectedZ.fromS.R(this.S, this.SST, this.imap, win, rratio, MIN.SNPs)
    if(is.null(this.Z)){
       return(list(bestchpt=c(NA,NA), bestZ = NA, Z=NA))
        
    }
    
    this.Z = (this.Z - 1)/sqrt(2)    
 
    maxind = matrix.max(this.Z)
    bestchpt= c(maxind[1],maxind[1]+maxind[2])
    bestZ = this.Z[maxind[1],maxind[2]]
    
    list(bestchpt=bestchpt,bestZ=bestZ, Z=this.Z)
}

