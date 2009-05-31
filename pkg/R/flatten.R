`flatten` <-
function(z){
  n=prod(dim(z))
  z.flat = matrix(data=z,nrow=n, ncol=1)
}

