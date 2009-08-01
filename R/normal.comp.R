`normal.comp` <-
function(y,mu,sigma){
  deter = det(sigma)
  inv = solve(sigma)
  a = y-mu
  b = a%*%inv%*%a
  result = deter^(-1/2)*exp(-0.5*b)
  
  # to avoid 0
  if (result==0){
    result = 1e-20
  }
  result
}

