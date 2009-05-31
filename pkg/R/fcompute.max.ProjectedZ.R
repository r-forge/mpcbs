`fcompute.max.ProjectedZ` <-
function(this.S,this.SST,this.imap,win,rratio,MIN.SNPs, f=NULL){
    K=ncol(this.S) # number of platforms.
    this.T = nrow(this.S)
    
    if(is.na(rratio) || (length(rratio) != K)){
        rratio = rep(1,K)
    }
    
    win = min(win, this.T-1)
    
    temp = fscan.max(this.S, this.SST, this.imap=this.imap,MIN.SNPs=MIN.SNPs,rratio=rratio,f=f, use.Project.statistic=TRUE, verbose=FALSE)
    bestchpt = temp$seg
    bestZ = temp$maxZ
    
    list(bestchpt=bestchpt,bestZ=bestZ)
}

