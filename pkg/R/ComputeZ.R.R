`ComputeZ.R` <-
function(Y, T, win, ALPHA){
    dfnum <- 1;
    dfden <- T-2;



    log.alpha <- log(ALPHA)
    
    g <- function(u){
        i2 <- (u < 10 + log.alpha)
        r1 <- u
        r1[i2] <- r1[i2] * exp(u[i2]/2)/(ALPHA+exp(u[i2]/2))
        return( r1 )
    }
    

    N <- dim(y)[1]

    U=matrix(nrow=T, ncol=win, data=0);
    Z=U;


    for(i in 1:N){

        
        S=cumsum(y[i,]);
        SST = sum((y[i,]-S[T]/T)^2);
        
        for(k in 1:win){
            diff1 <- S[(k+1):T]-S[1:(T-k)]
                SSb = k*(diff1/k-S[T]/T)^2;
                SSb = SSb + (T-k)*((S[T]-diff1)/(T-k)-S[T]/T)^2;
                SSw = SST-SSb;
                U[1:(T-k),k] = (SSb/dfnum)/(SSw/dfden);
        }
        
        Z <- Z+g(U);  # Z[t,k]: change starting at t+1, ending at t+k.
    }
    return(Z)
}

