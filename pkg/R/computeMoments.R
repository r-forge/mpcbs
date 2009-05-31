`computeMoments` <-
function(g,theta){

    INTLIM.THRESH=0.001
    psidottop.int <-function(u){ g(u)*exp(theta*g(u))*exp(-(u)^2/2) }    
    for( INTLIM in 10:50 ){
        temp=psidottop.int(INTLIM)
        if (is.na(temp)){
            INTLIM=INTLIM-1
            break
        }
        if (temp<INTLIM.THRESH) break
    }
    psidottop =  integrate(psidottop.int,lower=-INTLIM,upper=INTLIM)
    psidottop = psidottop$value/sqrt(2*pi)

    psidotbot.int =  function(u){ exp(theta*g(u))*exp(-(u)^2/2)}  
    for( INTLIM in 10:50 ){
        temp=psidotbot.int(INTLIM)
        if( is.na(temp)){
            INTLIM=INTLIM-1
            break
        }
        if( temp<INTLIM.THRESH ) break
    }
    psidotbot = integrate(psidotbot.int, lower=-INTLIM, upper=INTLIM)
    psidotbot = psidotbot$value/sqrt(2*pi)
    psidot = psidottop/psidotbot
    
    EgUsq.int<-function(u){ g(u)^2*exp(theta*g(u))*exp(-(u)^2/2) }    
    for( INTLIM in 10:50 ){
        temp=EgUsq.int(INTLIM)
        if( is.na(temp)){
            INTLIM=INTLIM-1
            break
        }
        if( temp<INTLIM.THRESH) break
    }
    EgUsq =  integrate(EgUsq.int,lower=-INTLIM,upper=INTLIM)
    EgUsq = EgUsq$value/sqrt(2*pi)

    psidotdot = (EgUsq*psidotbot - psidottop^2)/(psidotbot^2);
    psi = psidotbot

    list(psi=psi,psidot=psidot,psidotdot=psidotdot)
}

