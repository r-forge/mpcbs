`cluster.1` <-
function(mat1,mat2){
  mat = rbind(mat1,mat2)
  temp = getmusigma(mat)
  mu = temp$mu
  sigma = temp$sigma
  n = length(mat[,1])
  g.theta = 0
  for (i in 1:n){
    g.theta = g.theta + log(normal.comp(mat[i,],mu,sigma))
  }
  return(g.theta-3*log(n))
}

