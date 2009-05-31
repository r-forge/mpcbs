`ComputeBYZ.fromS.R.partial` <-
function(this.S,y.var,this.imap,start.inds, end.inds,delta,MIN.SNPs){
    T = nrow(this.S)
    N = ncol(this.S)
    dfnum <- 1;

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
          
          l.delta.1 = delta[1]*(diff1 - this.S[T,]*nsnps/totalsnps)/sqrt(y.var) - (delta[1]^2/2)*nsnps*(1-nsnps/totalsnps)
          l.delta.2 = delta[2]*(diff1 - this.S[T,]*nsnps/totalsnps)/sqrt(y.var) - (delta[2]^2/2)*nsnps*(1-nsnps/totalsnps)
          
          set.to.zero = which(nsnps<MIN.SNPs | ((totalsnps-nsnps)<MIN.SNPs))
          l.delta.1[set.to.zero] = 0
          l.delta.2[set.to.zero] = 0
          Z[i,j] = sum(pmax(l.delta.1,0)+pmax(l.delta.2,0))
        }
      }
    }

    Z
}

