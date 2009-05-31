`computeZ.onechange` <-
function(y,t,y.var){
  T =ncol(y)
# cat("In computeZ.onechange: T=", T," t=", t," nrow(y)=", nrow(y),"\n")

  St = apply(y[,1:t],1,sum)
  ST = apply(y,1,sum)
  Zsq = (St - (t/T)*ST)^2/(y.var*t*(1-t/T))
  sumZ=sum(Zsq)
  pval = 1-pchisq(sumZ,nrow(y))
#  cat("In computeZ.onechange: T=", T," t=", t," nrow(y)=", nrow(y)," sumZ=", sumZ,", pval=", pval,"\n")
  # Looks like the Z-values are immense, and rarely ever accepts the null.
#  cat(Zsq)
#  cat("\n\n")
  
  list(sumZ=sumZ, pval=pval)
}

