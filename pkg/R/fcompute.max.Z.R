`fcompute.max.Z` <-
function(this.y,win,y.var,ALPHA,MIN.SNPs,SINGLECHANGE.THRESH=0.0001){
    N=ncol(this.y)
    this.T=nrow(this.y)
    if(this.T<3*MIN.SNPs){
        
        this.Z = NA
        bestZ=NA
        bestchpt=c(NA,NA)

    } else {

        win = min(win, this.T-1)

        this.S = apply(this.y,2,cumsum)
        this.SST = apply(this.y,2,var)*(this.T-1)
        
        # Find [t1, t2] that maximizes Z using filtered scan.
        temp = fscan.max(this.S,this.SST,this.imap=NULL,MIN.SNPs=MIN.SNPs,ALPHA=ALPHA)
        bestchpt= temp$seg
        bestZ = temp$maxZ
        
        # Test the left change-point individually.
        if(bestZ!=0){
        
            if(bestchpt[1]==1){ 
                pval.L = 1
            } else {
                pval.L<-computeZ.onechange(t(this.y[1:bestchpt[2],]),bestchpt[1],y.var)$pval
            }
            if(bestchpt[2]==this.T){
                pval.R =1
            } else {
                pval.R<-computeZ.onechange(t(this.y[(bestchpt[1]+1):this.T,]),bestchpt[2]-bestchpt[1],y.var)$pval
            }
         # Pruning of left and right change-points (Olshen and Venkatraman suggestion) happens here, even though we haven't yet tested whether the double change-point has significant p-value.  The point is, in case the double change-point has significant p-value, then we would need to test for the significance of the left and right change-points anyway.
          if(pval.L>SINGLECHANGE.THRESH || bestchpt[1]<2*MIN.SNPs){
            # Added 6/15: we prune out change-points that are within MIN.SNPs to either end point.  In future, need
            # get rid of ##*## above and just use this step to control the boundary effects.
            
             #cat("Pruned out left change-point ", bestchpt[1]," out of ", this.T,"\n")
            bestchpt = c(bestchpt[2], NA)
          } else {
            if(pval.R>SINGLECHANGE.THRESH || (this.T-bestchpt[2]) < 2*MIN.SNPs){
#             cat("Pruned out right change-point ", bestchpt[2]," out of ", this.T,"\n")
             bestchpt = c(bestchpt[1], NA)
            }
          }
        } else {
          bestchpt = c(NA,NA)
        }
    }
    list(bestchpt=bestchpt,bestZ=bestZ)
}

