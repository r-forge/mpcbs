`pmarg.sumweightedchisq` <-
function(b,ALPHA,N){
    g<-function(u){
         u^2*exp(u^2/2)/(ALPHA+exp(u^2/2))
    }
    g.moments=computeMoments(g,0);
    gnormed<-function(u){(g(u)-g.moments$psidot)/sqrt(g.moments$psidotdot)}
    theta0=sqrt(g.moments$psidotdot)/2   # need theta/sqrt(psidotdot) < 1/2 for psi(theta)<infty.
    b0=b/sqrt(N)
    THRESH=0.0001
    marg = computeTiltDirect(b0,gnormed,THRESH,theta0)
    thetabN = marg$theta*sqrt(N)
    pmarg = exp(-thetabN*b + N*log(marg$psi))/sqrt(2*pi*marg$psidotdot)
    pmarg
}

