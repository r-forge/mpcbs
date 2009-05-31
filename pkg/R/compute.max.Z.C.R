`compute.max.Z.C` <-
function(this.y,win,y.var,ALPHA,MIN.SNPs,SINGLECHANGE.THRESH=0.0001){
    N=ncol(this.y)
    this.T=nrow(this.y)
    if(this.T<3*MIN.SNPs){
        
        this.Z = NA
        bestZ=NA
        bestchpt=c(NA,NA)

    } else {

        win = min(win, this.T-1)
        
        this.Z<-ComputeZ.C(t(this.y), this.T, win, ALPHA)    
        g2<-function(u){ u^2*exp(u^2/2)/(ALPHA+exp(u^2/2))}
        g2.moments=computeMoments(g2,0);
        this.Z = (this.Z - N*g2.moments$psidot)/sqrt(g2.moments$psidotdot*N)    
    

        ##*##
        this.Z[,1:(MIN.SNPs-1)] = 0
        this.Z[1:(MIN.SNPs-1),]=0
        for(i in 1:this.T-MIN.SNPs-1){
            maxwin = this.T-MIN.SNPs-i
            if(maxwin<win) this.Z[i, maxwin:win]=0
        }
        this.Z[(this.T-MIN.SNPs):this.T,]=0
        ##*##
        

        maxind = matrix.max(this.Z)
        bestchpt= c(maxind[1]+1,maxind[1]+maxind[2]+1)
        bestZ = this.Z[maxind[1],maxind[2]] # keep bestZ the best of the two-change-point scan.  If pruning out the left or right change-point made a big difference in bestZ, it would not have been pruned out anyway.    
        
        # Test the left change-point individually.
        if(bestZ!=0){
          pval.L<-computeZ.onechange(t(this.y[1:bestchpt[2],]),bestchpt[1],y.var)$pval
          pval.R<-computeZ.onechange(t(this.y[(bestchpt[1]+1):this.T,]),bestchpt[2]-bestchpt[1],y.var)$pval

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
    list(bestchpt=bestchpt,bestZ=bestZ, Z=this.Z)
}

