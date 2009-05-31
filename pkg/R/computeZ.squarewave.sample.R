`computeZ.squarewave.sample` <-
function(this.y,seg,y.var){
  T =ncol(this.y)
  Ss = apply(as.matrix(this.y[,1:seg[1]],nrow=nrow(this.y),ncol=seg[1]),1,sum)
  St = apply(this.y[,1:seg[2]],1,sum)
  ST = apply(this.y,1,sum)
  k=seg[2]-seg[1]
  Zsq = (St-Ss - (k/T)*ST)^2/(y.var*k*(1-k/T))
  pval = 1-pchisq(Zsq,1)

  list(Zsq=Zsq, pval=pval)
}

