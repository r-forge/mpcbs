`ComputeBYZ.fromS.R` <-
function(this.S,y.var,this.imap,win,delta,MIN.SNPs){
    T = nrow(this.S) # number of snps.
    N = ncol(this.S) # number of samples.
    dfnum <- 1;

    win = min(win, T-1)
    if(win==0){ 
        return(NULL)
    }

    l.delta.1=matrix(nrow=T, ncol=win, data=0);
    l.delta.2=matrix(nrow=T, ncol=win, data=0);
    Z=l.delta.1;
    for(i in 1:N){ 
        totalsnps = this.imap[T,i]-this.imap[1,i]
        dfden = totalsnps-2
        for(k in 1:win){
                nsnps = this.imap[(k+1):T,i]-this.imap[1:(T-k),i]
                diff1 <- this.S[(k+1):T, i]-this.S[1:(T-k),i]
                
                l.delta.1[1:(T-k),k]  = delta[1]*(diff1 - this.S[T,i]*nsnps/totalsnps)/sqrt(y.var[i]) - (delta[1]^2/2)*nsnps*(1-nsnps/totalsnps)
                l.delta.2[1:(T-k),k]  = delta[2]*(diff1 - this.S[T,i]*nsnps/totalsnps)/sqrt(y.var[i]) - (delta[2]^2/2)*nsnps*(1-nsnps/totalsnps)
                              
                set.to.zero = which(nsnps<MIN.SNPs | ((totalsnps-nsnps)<MIN.SNPs))
                l.delta.1[set.to.zero,k] = 0
                l.delta.2[set.to.zero,k] = 0
        }
        Z <- Z+pmax(l.delta.1,0)+pmax(l.delta.2,0)
    }
    return(Z)
}

