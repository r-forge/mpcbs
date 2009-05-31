`Chisq.Contrib.C` <-
function(y, T, chpts){
    res1 <- .Call("ChisqContrib", as.vector(y), T, dim(y), as.integer(chpts[,1]), as.integer(chpts[,2])  )
    return(matrix(data=res1, nrow=nrow(chpts))  )
}

