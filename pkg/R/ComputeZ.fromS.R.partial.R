`ComputeZ.fromS.R.partial` <-
function(this.S,this.SST,this.imap,start.inds, end.inds,ALPHA,MIN.SNPs){
    T = nrow(this.S)
    N = ncol(this.S)
    dfnum <- 1;

    log.alpha <- log(ALPHA)
    g <- function(u){
        i2 <- (u < 10 + log.alpha)
        r1 <- u
        r1[i2] <- r1[i2] * exp(u[i2]/2)/(ALPHA+exp(u[i2]/2))
        return( r1 )
    }
    
    Z = matrix(nrow=length(start.inds),ncol=length(end.inds),data=0)

    totalsnps = this.imap[T,] - this.imap[1,]
    dfden = totalsnps-2
    for(i in 1:length(start.inds)){
      for(j in 1:length(end.inds)){
        st = start.inds[i]
        ed = end.inds[j]
        if(st > 0 && st < ed && ed<T){
          nsnps = this.imap[ed,] - this.imap[st,]
          diff1 = this.S[ed,]-this.S[st,]
          SSb = nsnps*(diff1/nsnps - this.S[T,]/totalsnps)^2
          SSb = SSb + (totalsnps-nsnps)*((this.S[T,]-diff1)/(totalsnps-nsnps)-this.S[T,]/totalsnps)^2;
          SSw = this.SST - SSb
          U = (SSb/dfnum)/(SSw/dfden)
          set.to.zero = which(nsnps<MIN.SNPs | ((totalsnps-nsnps)<MIN.SNPs))
          U[set.to.zero] = 0
          Z[i,j] = sum(g(U))
        }
      }
    }

    Z
}

