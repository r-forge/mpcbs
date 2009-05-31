`computeZ.onechange.sample` <-
function(y,t,y.var){

  T =ncol(y)
#  cat("here!! T=", T,", t=", t,", nrow(y)=", nrow(y),"\n\n")

  St = apply(y[,1:t],1,sum)
  ST = apply(y,1,sum)
  Zsq = (St - (t/T)*ST)^2/(y.var*t*(1-t/T))
  pval = 1-pchisq(Zsq,1)

  list(Zsq, pval=pval)
}

