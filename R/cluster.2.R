`cluster.2` <-
function(mat1,mat2){
  temp1 = getmusigma(mat1)
  mu1 = temp1$mu
  sigma1 = temp1$sigma
  temp2 = getmusigma(mat2)
  mu2 = temp2$mu
  sigma2 = temp2$sigma
  n1 = length(mat1[,1])
  n2 = length(mat2[,1])
  g.theta = 0
  for (i in 1:n1){
    temp1 = normal.comp(mat1[i,],mu1,sigma1)
    temp2 = normal.comp(mat1[i,],mu2,sigma2)
    g.theta = g.theta + log(temp1+temp2)
  }
  for (i in 1:n2){
    temp1 = normal.comp(mat2[i,],mu1,sigma1)
    temp2 = normal.comp(mat2[i,],mu2,sigma2)
    g.theta = g.theta + log(temp1+temp2)
  }
  return(g.theta-(n1+n2)*log(2)-6*log(n1+n2))
}

