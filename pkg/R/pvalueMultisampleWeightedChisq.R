`pvalueMultisampleWeightedChisq` <-
function(b,m,delta,win,ALPHA,N){
    delta1 = min(1, win/m);
#    if(msscan.debug.trace){
#        print(paste("m = ", m, "; b = ", b, "; delta = ", delta, "; delta1 = ", delta1))
#    }

    beta=computeBeta(ALPHA)

    integrand<-function(u){
        vu(sqrt(2)*b*sqrt(beta)/(sqrt(m)*sqrt(u*(1-u))))^2/(u^2*(1-u))
    }
    
    integral =  integrate(integrand,lower=delta,upper=delta1)
    pmarg = pmarg.sumweightedchisq(b,ALPHA,N)

    pval=b^3*beta^2*pmarg*integral$value
    pval
}

