`computeTiltDirect` <-
function(b,g,THRESH,theta0) {

    theta = theta0 # if start theta at 0, then dh(theta)=0 at the first step for g(u)=u^2.
    prevtheta=Inf
    prevprevtheta=Inf
    thetarec = theta
    
    while( abs(theta-prevtheta)>THRESH && abs(theta-prevprevtheta)>THRESH ){
        
        g.moments<-computeMoments(g,theta)
        htheta=g.moments$psidot-b
        dhtheta = g.moments$psidotdot

        # cat("theta=", theta," htheta=",htheta,".\n",sep="")
        thetarec= c(thetarec, theta)
        prevprevtheta=prevtheta
        prevtheta=theta
        theta = prevtheta - htheta/dhtheta
    }

    theta=prevtheta
    psi = g.moments$psi

    list(theta=theta,psi=psi,psidot=g.moments$psidot,psidotdot=g.moments$psidotdot)
}

