`getmusigma` <-
function(mat){
  n = dim(mat)[1]
  mu = rep(0,2)
  mu[1] = mean(mat[,1])
  mu[2] = mean(mat[,2])
  sigma = matrix(0,2,2)
  for (i in 1:n){
    a = mat[i,]-mu
    sigma = sigma + a%*%t(a)
  }
  sigma = sigma/n
  list(mu=mu,sigma=sigma)
}

