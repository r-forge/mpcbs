`ComputeProjectedZ.fromS.R.partial` <-
function(this.S,this.SST,this.imap,start.inds, end.inds,rratio,MIN.SNPs){
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
          temp = diff1/nsnps - this.S[T,]/totalsnps
          SSb = nsnps*(temp)^2
          SSb = SSb + (totalsnps-nsnps)*((this.S[T,]-diff1)/(totalsnps-nsnps)-this.S[T,]/totalsnps)^2;
          SSw = this.SST - SSb
          U = sqrt((SSb/dfnum)/(SSw/dfden))
          
          weightU = sign(temp)*sqrt(nsnps*(1-nsnps/totalsnps)/(SSw/dfden))
#          weightU = sign(temp)
          
          set.to.zero = which(nsnps<MIN.SNPs | ((totalsnps-nsnps)<MIN.SNPs))
          U[set.to.zero] = 0
          Z[i,j] = sum(U*weightU*rratio)/sqrt(sum(rratio^2*weightU^2))
        }
      }
    }
    
    Z<- Z^2
    Z[which(is.na(Z), arr.ind=TRUE)]=0
 
    Z
}

