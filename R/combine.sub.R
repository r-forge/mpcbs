`combine.sub` <-
function(cnvobj, alpha=0.05){
  test = ABtest(cnvobj, alpha=alpha)
  if (!is.null(cnvobj$chpts.new)){
    chpts = cnvobj$chpts.new
    n = length(chpts)    
  }else{
    chpts = cnvobj$chpts
    n = length(chpts)
  }
  keep = rep(1,n)
  if (n>0){
    N = length(cnvobj$intensity[,1])
    chptall = c(1,chpts,(N+1))
    for (k in 1:n){
      if (test$major.result[2,k]==0 && test$minor.result[2,k]==0 && test$total.result[2,k]==0){
        keep[k] = 0
      }
    }
  }
  cnvobj$chpts.new = chpts[which(keep==1)]
  cnvobj
}

