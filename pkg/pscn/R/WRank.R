`WRank` <-
function(vec1, vec2, alpha=0.05){
  n = length(vec1)
  m = length(vec2)
  vec = c(vec1,vec2)
  Ws = sum(order(order(vec))[1:n])
  qvalue = (Ws-as.double(n)*(n+m+1)/2)/sqrt(as.double(m)*n*(n+m+1)/12)
  pvalue = pnorm(qvalue)
  result = 0
  if (pvalue<(alpha/2)){result = -1}
  if (pvalue>(1-alpha/2)){result = 1}
  return(c(pvalue,result))
}

