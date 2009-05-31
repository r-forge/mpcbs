`computeZ.test` <-
function(){
    par(mfrow=c(1,3))


    y <- matrix(rnorm(100*20), nrow=20)

    T = 100;
    win = 20;

 #   system.time(z <- ComputeZ.R(y, T, win, 0))
    system.time(z2 <- ComputeZ.C(y, T, win, 0))
    
    print(range(abs(z-z2))  )

    image(z)
    image(z2)
    image(z-z2)
    
}

