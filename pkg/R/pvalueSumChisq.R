`pvalueSumChisq` <-
function(b.normed,T,delta,T0,N){

    b=sqrt(b.normed*sqrt(2*N)+N)
    
    delta1 = min(1, T0/T)
    CbN = 1-(N-1)/b^2
    
    integrand<-function(u){
        vu(b*CbN/(sqrt(T)*sqrt(u*(1-u))))^2/(u^2*(1-u))
    }
    
    integral =  integrate(integrand,lower=delta,upper=delta1)
    
    pval=b^3*dchi(b,N)*CbN^3*2^{-2}*integral$value
    pval
}

