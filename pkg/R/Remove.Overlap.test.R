`Remove.Overlap.test` <-
function(){
    par(mfrow=c(1,3))

    y <- matrix(runif(2*20), ncol=2)
    y[,2] <- rowSums(y)

    i2a <- Remove.Overlap.R(y)
    i2b <- Remove.Overlap.C(y)

    

    print( sum(i2a != i2b) )

    i2a <- matrix(as.integer(i2a), nrow=1)
    i2b <- matrix(as.integer(i2a), nrow=1)

    image(i2a)
    image(i2b)
    image(i2a - i2b)
    
}

