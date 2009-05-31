`ComputeProjectedZ.fromS.R` <-
function(this.S,this.SST,this.imap,win,rratio,MIN.SNPs){
    T = nrow(this.S) # Number of SNPs.
    N = ncol(this.S) # Number of samples/platforms

    win = min(win, T-1)
    if(win==0){ 
        return(NULL)
    }
    
    U=matrix(nrow=T, ncol=win, data=0);
    weightU = U;
    sum.weight.sq = 0
    Z=U;
    for(i in 1:N){
        totalsnps = this.imap[T,i]-this.imap[1,i]
        dfden = totalsnps-2
        for(k in 1:win){
                nsnps = this.imap[(k+1):T,i]-this.imap[1:(T-k),i]
                diff1 <- this.S[(k+1):T, i]-this.S[1:(T-k),i]
                
                temp<-diff1/nsnps-this.S[T,i]/totalsnps
                
                SSb = nsnps*(temp)^2;
                SSb = SSb + (totalsnps-nsnps)*((this.S[T,i]-diff1)/(totalsnps-nsnps)-this.S[T,i]/totalsnps)^2;
                SSw = this.SST[i]-SSb;
                U[1:(T-k),k] = sqrt(SSb/(SSw/dfden));
                weightU[1:(T-k),k] = sign(temp)*sqrt(nsnps*(1-nsnps/totalsnps)/(SSw/dfden))
#                weightU[1:(T-k),k] = sign(temp)
                
                set.to.zero = which(nsnps<MIN.SNPs | ((totalsnps-nsnps)<MIN.SNPs))
                U[set.to.zero,k] = 0
        }
        Z <- Z+rratio[i]*weightU*U;
        sum.weight.sq = sum.weight.sq + rratio[i]^2*weightU^2
    }
    Z<- Z^2/sum.weight.sq
    Z[which(is.na(Z), arr.ind=TRUE)]=0
    return(Z)
}

