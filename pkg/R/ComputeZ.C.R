`ComputeZ.C` <-
function(y, T, win, ALPHA){
    res1 <- .Call("ComputeZ", as.vector(y), T, win, dim(y), ALPHA);

    Z <- matrix(nrow=T, ncol=win, data=res1)
    return(Z)
}

