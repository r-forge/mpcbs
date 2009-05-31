`fscan` <-
function(y,b,f=NULL,rho=0.99,MIN.SNPS=1, ALPHA=0, delta=c(1,-1), use.BY.statistic=FALSE, throw.small.Z=FALSE, MIN.ABS.MEDIAN=0.2, plots.pdf="fscanplots.pdf", verbose=TRUE){
  T=nrow(y)   # Number of SNPs
  N=ncol(y)   # Number of samples

  cat("Starting fscan with ",T," SNPs, ",N, " samples, threshold of ",b," and f=",f,"\n")

  if(is.null(f)){
  # This part doesn't work very well yet, initially, set f manually, e.g. f=c(0.01,0.1,1) works well.
    R= ceiling(log(log(T)) - log(-log(rho)))
    f = T^(-0.5^c(1:(R-1)))
    truncated.len = length(f)
    for(i in 2:length(f)){
      if(floor(1/f[i]) == floor(1/f[i-1])){
        truncated.len = i-1
        break
      }
    }
    f = c(f[1:truncated.len],1)
  }

  f = sort(f)
  if(f[length(f)] != 1) f = c(f,1)
  
  R = length(f)
  L = ceiling(f[2:R]/f[1:(R-1)])
  L = c(ceiling(T*f[1]),L)

  S = apply(y,2,cumsum)
  SST = apply(y,2,var)*(T-1)
  y.var <- compute.var(y)
 
  chpts = matrix(ncol=2,nrow=0)
  chpts.Z = matrix(nrow=1,ncol=0)

  g2<-function(u){ u^2*exp(u^2/2)/(ALPHA+exp(u^2/2))}
  g2.moments=computeMoments(g2,0);

  pdf(plots.pdf)
  
  for(r in 1:R){
    stepsize = floor(1/f[r])
    t = seq(1,T,stepsize) # t is the filtered anchor set.
    if(t[length(t)]<T) t = c(t,T) # always include the last datapoint in the set.

    this.S = S[t,]
    this.imap = matrix(rep(t,N),ncol=N,byrow=FALSE)

    # produce a diagnostic plot.
    plot(S[,1],xlab="SNP #", ylab="Cumulative sum of sample 1", main=paste("Round ",r,": Take 1 out of every ", stepsize," points.",sep=""))
    lines(t,this.S[,1],col="red")
    points(t,this.S[,1],col="red",pch=17,cex=2)

    # Refine previously found change-points using the denser anchor set.
    if(nrow(chpts)>0){
      for(i in 1:nrow(chpts)){
        ind.L = chpts[i,1] %/% stepsize
        ind.R = chpts[i,2] %/% stepsize
        
        check.win = f[r]/f[r-1]
        start.inds = c((ind.L - check.win):(ind.L+check.win))  
        end.inds = c((ind.R-check.win):(ind.R+check.win))
        
        if(use.BY.statistic){
            Z.part=ComputeBYZ.fromS.R.partial(this.S,y.var,this.imap,start.inds, end.inds,delta,MIN.SNPS)
            Z.part = Z.part/N
        } else {
            Z.part=ComputeZ.fromS.R.partial(this.S,SST,this.imap,start.inds, end.inds,ALPHA,MIN.SNPS)
            Z.part = (Z.part - N*g2.moments$psidot)/sqrt(g2.moments$psidotdot*N)          
        }

        #  image.plot(start.inds,end.inds,Z.part, xlab="Start of change - 1", ylab="End of change")

        maxind = matrix.max(Z.part)
        improved.cp = c(t[start.inds[maxind[1]]], t[end.inds[maxind[2]]])
        improved.Z = Z.part[maxind[1],maxind[2]]
        if(verbose) cat("Changepoints (",chpts[i,1],", ",chpts[i,2],") refined to (", improved.cp[1],", ",improved.cp[2],").  Z-score improved from ", chpts.Z[i]," to ",improved.Z,"\n", sep="")
        chpts[i,] = improved.cp
        chpts.Z[i] = improved.Z
      }
    }


    # Do scan with min window size MIN.SNPS and max window size L_i,    
    if(use.BY.statistic){
        Z = ComputeBYZ.fromS.R(this.S, y.var, this.imap, L[r], delta, MIN.SNPS)
        Z = Z/N
    } else {
        Z = ComputeZ.fromS.R(this.S, SST, this.imap, L[r], ALPHA, MIN.SNPS)
        Z = (Z - N*g2.moments$psidot)/sqrt(g2.moments$psidotdot*N)    
    }


    # produce a diagnostic plot.
    image.plot(t,c(1:L[r]),Z,xlab="Start of change - 1", ylab="Window size", main=paste("Sum of chisquare (Z-score) for round ",r,sep=""))
    
    Z.passinds = which(Z>b,arr.ind=TRUE)

    if(length(Z.passinds) > 0){

      # Some locations passed the threshold!  Gather together all that passed,
      # remove the (abundant) overlaps, and merge with changepoints found in previous rounds.
      
      Z.pass = Z[Z.passinds]
      Z.pass.order =order(Z.pass,decreasing=TRUE)
      Z.pass = Z.pass[Z.pass.order]
      Z.passinds = Z.passinds[Z.pass.order,]
    
      newchpts = cbind(Z.passinds[,1],Z.passinds[,1]+Z.passinds[,2])
      if(length(newchpts) == 2) newchpts =matrix(nrow=1,ncol=2,data=newchpts)
      newchpts.keep = Remove.Overlap.C(newchpts)
      newchpts = newchpts[newchpts.keep,]
      newchpts = matrix(t[newchpts],ncol=2,byrow=FALSE)
      newZ = Z.pass[newchpts.keep]

      # This is specific to the scanning algorithm and is not necessary if CBS were used:
      # If, say, in a region 1-1000, the subregion 500-1000 has an amplification.  Then
      # without knowledge of a "baseline", 1-500 would look like a deletion (with very high
      # Z-score.  It would get called at this stage, and all subsequent focal changes within
      # 1-500 would be missed.  To guard against this, conduct a simple filter here.  If
      # none of the samples have a |median| > MIN.ABS.MEDIAN, then throw this change-point away.
      # In fact, if a large portion of the chromosome were changed, then one should not use the
      # simple scanning algorithm but instead use CBS.
      throw=rep(0,nrow(newchpts))
      for(i in 1:nrow(newchpts)){
        y.med = apply(y[newchpts[i,1]:newchpts[i,2],], 2,median)
        if(sum(abs(y.med)>MIN.ABS.MEDIAN)==0) throw[i] = 1
      }
      newchpts=newchpts[!throw,]
      if(length(newchpts)==2)newchpts=matrix(nrow=1,ncol=2,data=newchpts)
      newZ=newZ[!throw]
      ####
      
      
      if(verbose) {
        cat("Change-points found that are new to scan ",r,":\n")
        print(cbind(newchpts,newZ), digits=1)
      }

      chpts = rbind(chpts,newchpts)
      chpts.Z = c(chpts.Z,newZ)
      chpts.len = chpts[,2]-chpts[,1]

      # Check again the updated list for overlaps with CNVs found in previous rounds,
      # and if there are any, discard the one with
      # the smaller Z score (or the smaller length).  This can be sped up!!
      if(throw.small.Z){ 
        chpts.Z.order = order(chpts.Z, decreasing=TRUE)
      } else {
        chpts.Z.order = order(chpts.len, decreasing=TRUE)
      }
      chpts.Z = chpts.Z[chpts.Z.order]
      chpts = chpts[chpts.Z.order,]
      if(length(chpts) == 2) chpts =matrix(nrow=1,ncol=2,data=chpts)
      chpts.keep = Remove.Overlap.C(chpts)
      chpts = chpts[chpts.keep,]
      if(length(chpts) == 2) chpts =matrix(nrow=1,ncol=2,data=chpts)      
      chpts.Z = chpts.Z[chpts.keep]
      
    }

    if(verbose){
        cat("After merging, change-points found so far:\n")
        if(length(chpts.Z)>0) print(cbind(chpts, chpts.Z), digits=1)
        cat("\n")
    }
  }

  dev.off()
  
  list(chpts = chpts, chpts.Z = chpts.Z)
}

