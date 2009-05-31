`ComputeZ.fromS.R` <-
function(this.S,this.SST,this.imap,win,ALPHA,MIN.SNPs){
    T = nrow(this.S)
    N = ncol(this.S)
    dfnum <- 1;

    win = min(win, T-1)
    if(win==0){ 
        return(NULL)
    }

    log.alpha <- log(ALPHA)
    g <- function(u){
        i2 <- (u < 10 + log.alpha)
        r1 <- u
        r1[i2] <- r1[i2] * exp(u[i2]/2)/(ALPHA+exp(u[i2]/2))
        return( r1 )
    }
    
    U=matrix(nrow=T, ncol=win, data=0);
    Z=U;
    for(i in 1:N){
        totalsnps = this.imap[T,i]-this.imap[1,i]
        dfden = totalsnps-2
        for(k in 1:win){
                nsnps = this.imap[(k+1):T,i]-this.imap[1:(T-k),i]
                diff1 <- this.S[(k+1):T, i]-this.S[1:(T-k),i]
                
                SSb = nsnps*(diff1/nsnps-this.S[T,i]/totalsnps)^2;
                SSb = SSb + (totalsnps-nsnps)*((this.S[T,i]-diff1)/(totalsnps-nsnps)-this.S[T,i]/totalsnps)^2;
                SSw = this.SST[i]-SSb;
                U[1:(T-k),k] = (SSb/dfnum)/(SSw/dfden);
                
                set.to.zero = which(nsnps<MIN.SNPs | ((totalsnps-nsnps)<MIN.SNPs))
                U[set.to.zero,k] = 0
        }
        Z <- Z+g(U);
    }
    return(Z)
}

