`Chisq.Contrib.test` <-
function(){
    y <- matrix(runif(100*20), ncol=20)
    T <- 5

    chpts <- matrix(runif(2*20), ncol=2)
    chpts[,2] <- rowSums(chpts)
    chpts[,1] <- as.integer(chpts[,1] * 15) + 1
    chpts[,2] <- as.integer(chpts[,2] * 15) + 3

    chpts[chpts > ncol(y)] <- ncol(y)
    
    
    res1 <- Chisq.Contrib.R(y, T, chpts)
    res2 <- Chisq.Contrib.C(y, T, chpts)

    print( range(abs(res1-res2))  )

    par(mfrow=c(1,3))
    image(res1)
    image(res2)
    image(res1-res2)
}

