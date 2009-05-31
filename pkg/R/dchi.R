`dchi` <-
function(y,N){
    fy <-(1-N/2)*log(2) + (N-1)*log(y) - y^2/2 - lgamma(N/2)
    exp(fy)
}

