`Change.Points.R` <-
function(ix, T, lookup){

    from.row <- (ix[1:lookup] %% T)
    from.col <- as.integer( ix[1:lookup] / T) + 1
    
    v1 <- from.row + 1
    v2 <- from.col + from.row
    return(cbind(v1, v2))
}

`Chisq.Contrib.C` <-
function(y, T, chpts){
    res1 <- .Call("ChisqContrib", as.vector(y), T, dim(y), as.integer(chpts[,1]), as.integer(chpts[,2]), PACKAGE="mpcbs")
    return(matrix(data=res1, nrow=nrow(chpts))  )
}

`Chisq.Contrib.fromS.R` <-
function(S,SST,imap,chpts){
    T = nrow(S)
    N = ncol(S)
    dfnum <- 1;
    n.segs <- nrow(chpts)
    
    ret.m <- matrix(nrow=n.segs, ncol=N, data=0)
 
    totalsnps = imap[T,]-imap[1,]
    
    for(i in 1:n.segs){
        st <- chpts[i,1]
        ed <- chpts[i,2]
        w=ed-st
        
        for(j in 1:N){
            nsnps = imap[ed,j]-imap[st,j]        
            mn1 <- S[T,j]/totalsnps[j]
            SSb <- (S[ed,j]-S[st,j] -nsnps*mn1)^2 / (nsnps*(1-nsnps/totalsnps[j]));
            SSw = SST[j]-SSb;
            ret.m[i,j] = (SSb/dfnum)/(SSw/(totalsnps[j]-2));
        }

    }
    
    return(ret.m)
}

`Chisq.Contrib.R` <-
function(y, T, chpts){
    dfnum <- 1;
    dfden <- T-2;

    N <- nrow(y)
    ret.m <- matrix(nrow=nrow(chpts), ncol=N, data=0)
    # chisq holds the contribution of each sample to the aberration.

    n.segs <- nrow(chpts)

    for(i in 1:n.segs){
        st <- chpts[i,1];
        ed <- chpts[i,2];
        w <- ed-st; # Changed 11/4, previously w<-ed-st+1.

        for(j in 1:N){
            mn1 <- sum(y[j,])/T
            SST <- sum( (y[j,]- mn1)^2 );
            SSb <- (( sum(y[j,(st+1):ed])-w*mn1)^2) / (w*(1-w/T));
            SSw = SST-SSb;
            ret.m[i,j] = (SSb/dfnum)/(SSw/dfden);
        }

    }
    return(ret.m)
}

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

`compute.max.ProjectedZ` <-
function(this.S,this.SST,this.imap,win,rratio,MIN.SNPs){
    K=ncol(this.S) # number of platforms.
    
    if(is.na(rratio) || (length(rratio) != K)){
        rratio = rep(1,K)
    }
    
    
    this.Z = ComputeProjectedZ.fromS.R(this.S, this.SST, this.imap, win, rratio, MIN.SNPs)
    if(is.null(this.Z)){
       return(list(bestchpt=c(NA,NA), bestZ = NA, Z=NA))
        
    }
    
    this.Z = (this.Z - 1)/sqrt(2)    
 
    maxind = matrix.max(this.Z)
    bestchpt= c(maxind[1],maxind[1]+maxind[2])
    bestZ = this.Z[maxind[1],maxind[2]]
    
    list(bestchpt=bestchpt,bestZ=bestZ, Z=this.Z)
}


`compute.max.Z.C` <-
function(this.y,win,y.var,ALPHA,MIN.SNPs,SINGLECHANGE.THRESH=0.0001){
    N=ncol(this.y)
    this.T=nrow(this.y)
    if(this.T<3*MIN.SNPs){
        
        this.Z = NA
        bestZ=NA
        bestchpt=c(NA,NA)

    } else {

        win = min(win, this.T-1)
        
        this.Z<-ComputeZ.C(t(this.y), this.T, win, ALPHA)    
        g2<-function(u){ u^2*exp(u^2/2)/(ALPHA+exp(u^2/2))}
        g2.moments=computeMoments(g2,0);
        this.Z = (this.Z - N*g2.moments$psidot)/sqrt(g2.moments$psidotdot*N)    
    

        ##*##
        this.Z[,1:(MIN.SNPs-1)] = 0
        this.Z[1:(MIN.SNPs-1),]=0
        for(i in 1:this.T-MIN.SNPs-1){
            maxwin = this.T-MIN.SNPs-i
            if(maxwin<win) this.Z[i, maxwin:win]=0
        }
        this.Z[(this.T-MIN.SNPs):this.T,]=0
        ##*##
        

        maxind = matrix.max(this.Z)
        bestchpt= c(maxind[1]+1,maxind[1]+maxind[2]+1)
        bestZ = this.Z[maxind[1],maxind[2]] # keep bestZ the best of the two-change-point scan.  If pruning out the left or right change-point made a big difference in bestZ, it would not have been pruned out anyway.    
        
        # Test the left change-point individually.
        if(bestZ!=0){
          pval.L<-computeZ.onechange(t(this.y[1:bestchpt[2],]),bestchpt[1],y.var)$pval
          pval.R<-computeZ.onechange(t(this.y[(bestchpt[1]+1):this.T,]),bestchpt[2]-bestchpt[1],y.var)$pval

         # Pruning of left and right change-points (Olshen and Venkatraman suggestion) happens here, even though we haven't yet tested whether the double change-point has significant p-value.  The point is, in case the double change-point has significant p-value, then we would need to test for the significance of the left and right change-points anyway.
          if(pval.L>SINGLECHANGE.THRESH || bestchpt[1]<2*MIN.SNPs){
            # Added 6/15: we prune out change-points that are within MIN.SNPs to either end point.  In future, need
            # get rid of ##*## above and just use this step to control the boundary effects.
            
             #cat("Pruned out left change-point ", bestchpt[1]," out of ", this.T,"\n")
            bestchpt = c(bestchpt[2], NA)
          } else {
            if(pval.R>SINGLECHANGE.THRESH || (this.T-bestchpt[2]) < 2*MIN.SNPs){
#             cat("Pruned out right change-point ", bestchpt[2]," out of ", this.T,"\n")
             bestchpt = c(bestchpt[1], NA)
            }
          }
        } else {
          bestchpt = c(NA,NA)
        }
    }
    list(bestchpt=bestchpt,bestZ=bestZ, Z=this.Z)
}


`compute.max.Z` <-
function(this.S,this.SST,this.imap,win,ALPHA,MIN.SNPs){
    K=ncol(this.S)
    
    this.Z = ComputeZ.fromS.R(this.S, this.SST, this.imap, win, ALPHA, MIN.SNPs)
    if(is.null(this.Z)){
       return(list(bestchpt=c(NA,NA), bestZ = NA, Z=NA))
        
    }
    
    g2<-function(u){ u^2*exp(u^2/2)/(ALPHA+exp(u^2/2))}
    g2.moments=computeMoments(g2,0);
    this.Z = (this.Z - K*g2.moments$psidot)/sqrt(g2.moments$psidotdot*K)    
 
    maxind = matrix.max(this.Z)
    bestchpt= c(maxind[1],maxind[1]+maxind[2])
    bestZ = this.Z[maxind[1],maxind[2]]
    
    list(bestchpt=bestchpt,bestZ=bestZ, Z=this.Z)
 
}


`computeBeta` <-
function(ALPHA){
    w<-function(u){
         exp(u^2/2)/(ALPHA+exp(u^2/2))
    }
    integrand<-function(u){
        u^2*w(u)^2*(2*u^4*w(u)*(1-w(u))+5*u^2*w(u)-3*u^2-2)*dchi(u,1)
    }
    
#    us=seq(0,10,0.1)
#    ys=integrand(us)
#    plot(us,ys)
    
    numerator=integrate(integrand,lower=0,upper=10)
    
    g<-function(u){
         u^2*exp(u^2/2)/(ALPHA+exp(u^2/2))
    }
    g.moments=computeMoments(g,0)
    numerator$value/(2*g.moments$psidotdot)
}

`ComputeBYZ.fromS.R.partial` <-
function(this.S,y.var,this.imap,start.inds, end.inds,delta,MIN.SNPs){
    T = nrow(this.S)
    N = ncol(this.S)
    dfnum <- 1;

    Z = matrix(nrow=length(start.inds),ncol=length(end.inds),data=0)

    totalsnps = this.imap[T,] - this.imap[1,]
    dfden = totalsnps-2
    for(i in 1:length(start.inds)){
      for(j in 1:length(end.inds)){
        st = start.inds[i]
        ed = end.inds[j]
        if(st > 0 && st < ed && ed<T){
          nsnps = this.imap[ed,] - this.imap[st,]
          diff1 = this.S[ed,]-this.S[st,]
          
          l.delta.1 = delta[1]*(diff1 - this.S[T,]*nsnps/totalsnps)/sqrt(y.var) - (delta[1]^2/2)*nsnps*(1-nsnps/totalsnps)
          l.delta.2 = delta[2]*(diff1 - this.S[T,]*nsnps/totalsnps)/sqrt(y.var) - (delta[2]^2/2)*nsnps*(1-nsnps/totalsnps)
          
          set.to.zero = which(nsnps<MIN.SNPs | ((totalsnps-nsnps)<MIN.SNPs))
          l.delta.1[set.to.zero] = 0
          l.delta.2[set.to.zero] = 0
          Z[i,j] = sum(pmax(l.delta.1,0)+pmax(l.delta.2,0))
        }
      }
    }

    Z
}

`ComputeBYZ.fromS.R` <-
function(this.S,y.var,this.imap,win,delta,MIN.SNPs){
    T = nrow(this.S) # number of snps.
    N = ncol(this.S) # number of samples.
    dfnum <- 1;

    win = min(win, T-1)
    if(win==0){ 
        return(NULL)
    }

    l.delta.1=matrix(nrow=T, ncol=win, data=0);
    l.delta.2=matrix(nrow=T, ncol=win, data=0);
    Z=l.delta.1;
    for(i in 1:N){ 
        totalsnps = this.imap[T,i]-this.imap[1,i]
        dfden = totalsnps-2
        for(k in 1:win){
                nsnps = this.imap[(k+1):T,i]-this.imap[1:(T-k),i]
                diff1 <- this.S[(k+1):T, i]-this.S[1:(T-k),i]
                
                l.delta.1[1:(T-k),k]  = delta[1]*(diff1 - this.S[T,i]*nsnps/totalsnps)/sqrt(y.var[i]) - (delta[1]^2/2)*nsnps*(1-nsnps/totalsnps)
                l.delta.2[1:(T-k),k]  = delta[2]*(diff1 - this.S[T,i]*nsnps/totalsnps)/sqrt(y.var[i]) - (delta[2]^2/2)*nsnps*(1-nsnps/totalsnps)
                              
                set.to.zero = which(nsnps<MIN.SNPs | ((totalsnps-nsnps)<MIN.SNPs))
                l.delta.1[set.to.zero,k] = 0
                l.delta.2[set.to.zero,k] = 0
        }
        Z <- Z+pmax(l.delta.1,0)+pmax(l.delta.2,0)
    }
    return(Z)
}

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

`ComputeProjectedZ.fromS.R.segments`<-
function(this.S, this.SST,this.imap,segs,rratio,MIN.SNPs){
    T = nrow(this.S)
    N = ncol(this.S)
    dfnum <- 1;
    totalsnps = this.imap[T,] - this.imap[1,]
    dfden = totalsnps-2


    m = nrow(segs)

    Z = rep(0,m)
    X = matrix(nrow=m, ncol=length(rratio),data=0)
    w = matrix(nrow=m, ncol=length(rratio),data=0)

    for(i in 1:m){
        st = segs[i,1]
        ed = segs[i,2]
        if(st > 0 && st < ed && ed<T){

          nsnps = this.imap[ed,] - this.imap[st,]
          diff1 = this.S[ed,]-this.S[st,]
          temp = diff1/nsnps - this.S[T,]/totalsnps
          SSb = nsnps*(temp)^2
          SSb = SSb + (totalsnps-nsnps)*((this.S[T,]-diff1)/(totalsnps-nsnps)-this.S[T,]/totalsnps)^2;
          SSw = this.SST - SSb
          U = sqrt((SSb/dfnum)/(SSw/dfden))

          weightU = sign(temp)*sqrt(nsnps*(1-nsnps/totalsnps)/(SSw/dfden))

          weightU[which(is.na(weightU))]=0


          set.to.zero = which(nsnps<MIN.SNPs | ((totalsnps-nsnps)<MIN.SNPs))
          U[set.to.zero] = 0
          Z[i] = sum(U*weightU*rratio)/sqrt(sum(rratio^2*weightU^2))
          Z[i] = Z[i]^2
          X[i,] = U
          w[i,] = abs(weightU*rratio/sqrt(sum(rratio^2*weightU^2)))

        }
    }

    list(Z = Z, X = X, w = w)
}


`ComputeProjectedZ.fromS.R.partial` <-
function(this.S,this.SST,this.imap,start.inds, end.inds,rratio,MIN.SNPs){
    T = nrow(this.S)
    N = ncol(this.S)
    dfnum <- 1;


    Z = matrix(nrow=length(start.inds),ncol=length(end.inds),data=0)

    totalsnps = this.imap[T,] - this.imap[1,]
    dfden = totalsnps-2
    for(i in 1:length(start.inds)){
      for(j in 1:length(end.inds)){
        st = start.inds[i]
        ed = end.inds[j]
        if(st > 0 && st < ed && ed<T){
          nsnps = this.imap[ed,] - this.imap[st,]
          diff1 = this.S[ed,]-this.S[st,]
          temp = diff1/nsnps - this.S[T,]/totalsnps
          SSb = nsnps*(temp)^2
          SSb = SSb + (totalsnps-nsnps)*((this.S[T,]-diff1)/(totalsnps-nsnps)-this.S[T,]/totalsnps)^2;
          SSw = this.SST - SSb
          U = sqrt((SSb/dfnum)/(SSw/dfden))
          
          weightU = sign(temp)*sqrt(nsnps*(1-nsnps/totalsnps)/(SSw/dfden))
#          weightU = sign(temp)
          
          set.to.zero = which(nsnps<MIN.SNPs | ((totalsnps-nsnps)<MIN.SNPs))
          U[set.to.zero] = 0
          Z[i,j] = sum(U*weightU*rratio)/sqrt(sum(rratio^2*weightU^2))
        }
      }
    }
    
    Z<- Z^2
    Z[which(is.na(Z), arr.ind=TRUE)]=0
 
    Z
}


`ComputeProjectedZ.fromS.R` <-
function(this.S,this.SST,this.imap,win,rratio,MIN.SNPs){
    T = nrow(this.S) # Number of SNPs.
    N = ncol(this.S) # Number of samples/platforms

    win = min(win, T-1)
    if(win==0){ 
        return(NULL)
    }
    
    U=matrix(nrow=T, ncol=win, data=0);
    weightU = U;
    sum.weight.sq = 0
    Z=U;
    for(i in 1:N){
        totalsnps = this.imap[T,i]-this.imap[1,i]
        dfden = totalsnps-2
        for(k in 1:win){
                nsnps = this.imap[(k+1):T,i]-this.imap[1:(T-k),i]
                diff1 <- this.S[(k+1):T, i]-this.S[1:(T-k),i]
                
                temp<-diff1/nsnps-this.S[T,i]/totalsnps
                
                SSb = nsnps*(temp)^2;
                SSb = SSb + (totalsnps-nsnps)*((this.S[T,i]-diff1)/(totalsnps-nsnps)-this.S[T,i]/totalsnps)^2;
                SSw = this.SST[i]-SSb;
                U[1:(T-k),k] = sqrt(SSb/(SSw/dfden));
                weightU[1:(T-k),k] = sign(temp)*sqrt(nsnps*(1-nsnps/totalsnps)/(SSw/dfden))
#                weightU[1:(T-k),k] = sign(temp)
                
                set.to.zero = which(nsnps<MIN.SNPs | ((totalsnps-nsnps)<MIN.SNPs))
                U[set.to.zero,k] = 0
        }
        Z <- Z+rratio[i]*weightU*U;
        sum.weight.sq = sum.weight.sq + rratio[i]^2*weightU^2
    }
    Z<- Z^2/sum.weight.sq
    Z[which(is.na(Z), arr.ind=TRUE)]=0
    return(Z)
}

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


`ComputeZ.C` <-
function(y, T, win, ALPHA){
    res1 <- .Call("ComputeZ", as.vector(y), T, win, dim(y), ALPHA, PACKAGE="mpcbs");

    Z <- matrix(nrow=T, ncol=win, data=res1)
    return(Z)
}


`ComputeZ.fromS.R.partial` <-
function(this.S,this.SST,this.imap,start.inds, end.inds,ALPHA,MIN.SNPs){
    T = nrow(this.S)
    N = ncol(this.S)
    dfnum <- 1;

    log.alpha <- log(ALPHA)
    g <- function(u){
        i2 <- (u < 10 + log.alpha)
        r1 <- u
        r1[i2] <- r1[i2] * exp(u[i2]/2)/(ALPHA+exp(u[i2]/2))
        return( r1 )
    }
    
    Z = matrix(nrow=length(start.inds),ncol=length(end.inds),data=0)

    totalsnps = this.imap[T,] - this.imap[1,]
    dfden = totalsnps-2
    for(i in 1:length(start.inds)){
      for(j in 1:length(end.inds)){
        st = start.inds[i]
        ed = end.inds[j]
        if(st > 0 && st < ed && ed<T){
          nsnps = this.imap[ed,] - this.imap[st,]
          diff1 = this.S[ed,]-this.S[st,]
          SSb = nsnps*(diff1/nsnps - this.S[T,]/totalsnps)^2
          SSb = SSb + (totalsnps-nsnps)*((this.S[T,]-diff1)/(totalsnps-nsnps)-this.S[T,]/totalsnps)^2;
          SSw = this.SST - SSb
          U = (SSb/dfnum)/(SSw/dfden)
          set.to.zero = which(nsnps<MIN.SNPs | ((totalsnps-nsnps)<MIN.SNPs))
          U[set.to.zero] = 0
          Z[i,j] = sum(g(U))
        }
      }
    }

    Z
}


`ComputeZ.fromS.R` <-
function(this.S,this.SST,this.imap,win,ALPHA,MIN.SNPs){
    T = nrow(this.S)
    N = ncol(this.S)
    dfnum <- 1;

    win = min(win, T-1)
    if(win==0){ 
        return(NULL)
    }

    log.alpha <- log(ALPHA)
    g <- function(u){
        i2 <- (u < 10 + log.alpha)
        r1 <- u
        r1[i2] <- r1[i2] * exp(u[i2]/2)/(ALPHA+exp(u[i2]/2))
        return( r1 )
    }
    
    U=matrix(nrow=T, ncol=win, data=0);
    Z=U;
    for(i in 1:N){
        totalsnps = this.imap[T,i]-this.imap[1,i]
        dfden = totalsnps-2
        for(k in 1:win){
                nsnps = this.imap[(k+1):T,i]-this.imap[1:(T-k),i]
                diff1 <- this.S[(k+1):T, i]-this.S[1:(T-k),i]
                
                SSb = nsnps*(diff1/nsnps-this.S[T,i]/totalsnps)^2;
                SSb = SSb + (totalsnps-nsnps)*((this.S[T,i]-diff1)/(totalsnps-nsnps)-this.S[T,i]/totalsnps)^2;
                SSw = this.SST[i]-SSb;
                U[1:(T-k),k] = (SSb/dfnum)/(SSw/dfden);
                
                set.to.zero = which(nsnps<MIN.SNPs | ((totalsnps-nsnps)<MIN.SNPs))
                U[set.to.zero,k] = 0
        }
        Z <- Z+g(U);
    }
    return(Z)
}


`computeZ.onechange` <-
function(y,t,y.var){
  T =ncol(y)
# cat("In computeZ.onechange: T=", T," t=", t," nrow(y)=", nrow(y),"\n")

  St = apply(y[,1:t],1,sum)
  ST = apply(y,1,sum)
  Zsq = (St - (t/T)*ST)^2/(y.var*t*(1-t/T))
  sumZ=sum(Zsq)
  pval = 1-pchisq(sumZ,nrow(y))
#  cat("In computeZ.onechange: T=", T," t=", t," nrow(y)=", nrow(y)," sumZ=", sumZ,", pval=", pval,"\n")
  # Looks like the Z-values are immense, and rarely ever accepts the null.
#  cat(Zsq)
#  cat("\n\n")
  
  list(sumZ=sumZ, pval=pval)
}


`computeZ.onechange.sample` <-
function(y,t,y.var){

  T =ncol(y)
#  cat("here!! T=", T,", t=", t,", nrow(y)=", nrow(y),"\n\n")

  St = apply(y[,1:t],1,sum)
  ST = apply(y,1,sum)
  Zsq = (St - (t/T)*ST)^2/(y.var*t*(1-t/T))
  pval = 1-pchisq(Zsq,1)

  list(Zsq, pval=pval)
}


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

`computeZ.squarewave.sample` <-
function(this.y,seg,y.var){
  T =ncol(this.y)
  Ss = apply(as.matrix(this.y[,1:seg[1]],nrow=nrow(this.y),ncol=seg[1]),1,sum)
  St = apply(this.y[,1:seg[2]],1,sum)
  ST = apply(this.y,1,sum)
  k=seg[2]-seg[1]
  Zsq = (St-Ss - (k/T)*ST)^2/(y.var*k*(1-k/T))
  pval = 1-pchisq(Zsq,1)

  list(Zsq=Zsq, pval=pval)
}


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

`dchi` <-
function(y,N){
    fy <-(1-N/2)*log(2) + (N-1)*log(y) - y^2/2 - lgamma(N/2)
    exp(fy)
}

`fcompute.max.ProjectedZ` <-
function(this.S,this.SST,this.imap,win,rratio,MIN.SNPs, f=NULL){
    K=ncol(this.S) # number of platforms.
    this.T = nrow(this.S)
    
    if(is.na(rratio) || (length(rratio) != K)){
        rratio = rep(1,K)
    }
    
    win = min(win, this.T-1)
    
    temp = fscan.max(this.S, this.SST, this.imap=this.imap,MIN.SNPs=MIN.SNPs,rratio=rratio,f=f, use.Project.statistic=TRUE, verbose=FALSE)
    bestchpt = temp$seg
    bestZ = temp$maxZ
    
    list(bestchpt=bestchpt,bestZ=bestZ)
}


`fcompute.max.Z` <-
function(this.y,win,y.var,ALPHA,MIN.SNPs,SINGLECHANGE.THRESH=0.0001){
    N=ncol(this.y)
    this.T=nrow(this.y)
    if(this.T<3*MIN.SNPs){
        
        this.Z = NA
        bestZ=NA
        bestchpt=c(NA,NA)

    } else {

        win = min(win, this.T-1)

        this.S = apply(this.y,2,cumsum)
        this.SST = apply(this.y,2,var)*(this.T-1)
        
        # Find [t1, t2] that maximizes Z using filtered scan.
        temp = fscan.max(this.S,this.SST,this.imap=NULL,MIN.SNPs=MIN.SNPs,ALPHA=ALPHA)
        bestchpt= temp$seg
        bestZ = temp$maxZ
        
        # Test the left change-point individually.
        if(bestZ!=0){
        
            if(bestchpt[1]==1){ 
                pval.L = 1
            } else {
                pval.L<-computeZ.onechange(t(this.y[1:bestchpt[2],]),bestchpt[1],y.var)$pval
            }
            if(bestchpt[2]==this.T){
                pval.R =1
            } else {
                pval.R<-computeZ.onechange(t(this.y[(bestchpt[1]+1):this.T,]),bestchpt[2]-bestchpt[1],y.var)$pval
            }
         # Pruning of left and right change-points (Olshen and Venkatraman suggestion) happens here, even though we haven't yet tested whether the double change-point has significant p-value.  The point is, in case the double change-point has significant p-value, then we would need to test for the significance of the left and right change-points anyway.
          if(pval.L>SINGLECHANGE.THRESH || bestchpt[1]<2*MIN.SNPs){
            # Added 6/15: we prune out change-points that are within MIN.SNPs to either end point.  In future, need
            # get rid of ##*## above and just use this step to control the boundary effects.
            
             #cat("Pruned out left change-point ", bestchpt[1]," out of ", this.T,"\n")
            bestchpt = c(bestchpt[2], NA)
          } else {
            if(pval.R>SINGLECHANGE.THRESH || (this.T-bestchpt[2]) < 2*MIN.SNPs){
#             cat("Pruned out right change-point ", bestchpt[2]," out of ", this.T,"\n")
             bestchpt = c(bestchpt[1], NA)
            }
          }
        } else {
          bestchpt = c(NA,NA)
        }
    }
    list(bestchpt=bestchpt,bestZ=bestZ)
}

`flatten` <-
function(z){
  n=prod(dim(z))
  z.flat = matrix(data=z,nrow=n, ncol=1)
}

`fscan.max` <-
function(this.S, this.SST,  y.var=NULL, use.BY.statistic=FALSE, use.Project.statistic=FALSE, this.imap=NULL, 
                    f=NULL,MIN.SNPs=2,ALPHA=0,delta=c(1,-1), rratio=NULL, doplots=FALSE, verbose=FALSE){
  T=nrow(this.S)   # Number of SNPs
  N=ncol(this.S)   # Number of samples
  if(T<=MIN.SNPs){
    maxZ=NA
    seg=c(NA,NA)
  } else {

      if(is.null(f)){
        f.power = seq(min(-floor(log10(T))+1,0),0,1)
        f=10^f.power
      }
      
      if(T<1000) f=1 # even if a value for f is passed in, do not allow filtering for short sequences.
      
      f = sort(f)
      if(f[length(f)] != 1) f = c(f,1)
      R = length(f)
      L = ceiling(f[2:R]/f[1:(R-1)])
      L = c(ceiling(T*f[1]),L)
      
      chpts = matrix(ncol=2,nrow=0)
      chpts.Z = matrix(nrow=1,ncol=0)
    
      if(is.null(this.imap)){
        this.imap = matrix(rep(seq(1,T,1),N),ncol=N,byrow=FALSE)
      }
    
      g2<-function(u){ u^2*exp(u^2/2)/(ALPHA+exp(u^2/2))}
      g2.moments=computeMoments(g2,0);
      if(is.null(rratio)) rratio = rep(1,N)
      
      for(r in 1:R){
        stepsize = floor(1/f[r])
        
        t = seq(1,T,stepsize) # t is the filtered anchor set.
        if(t[length(t)]<T) t = c(t,T) # always include the last datapoint in the set.
    
        f.S = this.S[t,]
        f.imap = this.imap[t,]
    
        # produce a diagnostic plot.
        if(doplots){
          plot(this.S[,1],xlim=c(2000,3000),xlab="SNP #", ylab="Cumulative sum of sample 1", main=paste("Round ",r,": Take 1 out of every ", stepsize," points.",sep=""))
          lines(t,f.S[,1],col="red")
          points(t,f.S[,1],col="red",pch=17,cex=2)
        }
        
        # Refine previously found change-points using the denser anchor set.
        if(nrow(chpts)>0){
          for(i in 1:nrow(chpts)){
            ind.L = chpts[i,1] %/% stepsize
            ind.R = chpts[i,2] %/% stepsize
            
            check.win = f[r]/f[r-1]
	    start.inds = c((ind.L - check.win):(ind.L+check.win))  
            end.inds = c((ind.R-check.win):(ind.R+check.win))
            
                  
            if(use.BY.statistic){
                Z.part=ComputeBYZ.fromS.R.partial(f.S,y.var,f.imap,start.inds, end.inds,delta,MIN.SNPS)
                Z.part = Z.part/N
            } else {
                if(use.Project.statistic){
                    Z.part=ComputeProjectedZ.fromS.R.partial(f.S,this.SST,f.imap,start.inds, end.inds,rratio,MIN.SNPs)
                    Z.part = (Z.part -1)/sqrt(2)          
                } else {
                    Z.part=ComputeZ.fromS.R.partial(f.S,this.SST,f.imap,start.inds, end.inds,ALPHA,MIN.SNPs)
                    Z.part = (Z.part - N*g2.moments$psidot)/sqrt(g2.moments$psidotdot*N)          
                }
            }
    
            #  image.plot(start.inds,end.inds,Z.part, xlab="Start of change - 1", ylab="End of change")
    
            maxind = matrix.max(Z.part)
            improved.cp = c(t[start.inds[maxind[1]]], t[end.inds[maxind[2]]])
            improved.Z = Z.part[maxind[1],maxind[2]]
            if(verbose) cat("fscan.max: Changepoints (",chpts[i,1],", ",chpts[i,2],") refined to (", improved.cp[1],", ",improved.cp[2],").  Z-score improved from ", chpts.Z[i]," to ",improved.Z,"\n", sep="")
	    if(length(improved.cp)>=2){
	    chpts[i,] = improved.cp[1:2]
	    }else{
	    chpts[i,1]=improved.cp[1]
	    chpts[i,2]=improved.cp[1]
		}
            chpts.Z[i] = improved.Z
          }
        }
        
        
        # Do scan with min window size MIN.SNPS and max window size L_i,    
        if(use.BY.statistic){
            Z = ComputeBYZ.fromS.R(f.S, y.var, f.imap, L[r], delta, MIN.SNPS)
            Z = Z/N
        } else {
            if(use.Project.statistic){
                Z = ComputeProjectedZ.fromS.R(f.S, this.SST, f.imap, L[r], rratio, MIN.SNPs)
                Z = (Z - 1)/sqrt(2)    
            } else {
                # Do scan with min window size MIN.SNPs and max window size L_i,
                Z = ComputeZ.fromS.R(f.S, this.SST, f.imap, L[r], ALPHA, MIN.SNPs)
                Z = (Z - N*g2.moments$psidot)/sqrt(g2.moments$psidotdot*N)      
            }
        }
        # produce a diagnostic plot.
        if(doplots) image.plot(t,c(1:L[r]),Z, xlab="Start of change - 1", ylab="Window size", main=paste("Sum of chisquare (Z-score) for round ",r,sep=""))
    
        maxind = matrix.max(Z)
        #  newchpt= c(t[maxind[1]],t[maxind[1]+maxind[2]])
        newchpt= c(t[maxind[1]],t[maxind[1]+maxind[2]-1])  # changed 10/22, if maxind[1]+maxind[2] might become length(t)+1...
        newZ = Z[maxind[1],maxind[2]]
        if(verbose){
          cat("fscan.max: Change-point found in round ",r,":\n")
          print(c(newchpt,newZ), digits=1)
        }
        chpts = rbind(chpts,newchpt)
        if(length(chpts) == 2) chpts =matrix(nrow=1,ncol=2,data=chpts)      
        chpts.Z = c(chpts.Z,newZ)
      }
    
      # Take the maximum over the rounds, return chpt and Z score.
      max.ind = which.max(chpts.Z)
      maxZ = chpts.Z[max.ind]
      seg = chpts[max.ind,]
  }
  list(maxZ = maxZ, seg = seg)
}

`getCutoffEpidemicChangeLocationKnown` <-
function(pval){
    qchisq(1-pval,df=1)
}



`getCutoffMultisampleWeightedChisq` <-
function(pval,m,delta,win,N,alpha){
    
    cat("Computing threshold for weighted chi-square...\n")
    THRES = 0.1*pval
    currb = 1
    prevsmallerb = currb
    prevlargerb = 200
    currpval=1
    
    while( abs(currpval-pval)>THRES ){
#        cat("pval =", pval, ", currpval = ",currpval,", THRES=", THRES,".\n",sep="")
        
        if( currpval>pval){
            # need to increase b.
            prevsmallerb = currb
            currb = currb+(prevlargerb-currb)/2
        } else {
            # need to decrease b.
            prevlargerb = currb
            currb = currb - (currb-prevsmallerb)/2
        }
    
#        cat("currb = ",currb,"\n")
        currpval = pvalueMultisampleWeightedChisq(currb,m,delta,win,alpha,N);
    }    
    currb
}
`intmode` <-
function(x){
  xt<-tabulate(x)
  xmode<-which(xt == max(xt))
  xmode
}

`matrix.max` <-
function(Z){
    # If Z has NA values treat as -Inf
    na.inds = which(is.na(Z),arr.ind=TRUE)
    Z[na.inds]=-Inf

    row.max.col=max.col(Z, ties.method=c("random", "first", "last"))
    max.row = which.max(Z[cbind(1:nrow(Z),row.max.col)])
    c(max.row,row.max.col[max.row])
}



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

`pvalueMultisampleWeightedChisqold` <-
function(b,m,delta,win,ALPHA,N){
    delta1 = min(1, win/m);
#    if(msscan.debug.trace){
#        print(paste("m = ", m, "; b = ", b, "; delta = ", delta, "; delta1 = ", delta1))
#    }

    beta=computeBeta(ALPHA)

    integrand<-function(u){
        vu(sqrt(2)*b/(sqrt(m)*sqrt(u*(1-u))))^2/(u^2*(1-u))
    }
    
    integral =  integrate(integrand,lower=delta,upper=delta1)
    pmarg = pmarg.sumweightedchisq(b,ALPHA,N)

    pval=b^3*beta^2*pmarg*integral$value
    pval
}
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

`Remove.Overlap.C` <-
function(chpts){
    res1 <- .Call("RemoveOverlap", chpts[,1], chpts[,2], PACKAGE="mpcbs");
    
    i2 <- rep(FALSE, nrow(chpts))
    i2[res1 == 1] <- TRUE
    return(i2)
}







`Remove.Overlap.R` <-
function(chpts){

    lookup <- nrow(chpts)
    
    ret.v <- rep(FALSE, lookup)
    ret.v[1] <- TRUE
    for(i in 2:lookup){
        overlap <- 0;
        st <- chpts[i,1];
        ed <- chpts[i,2];
        for(j in which(ret.v)){
            # check for overlap.
            if(     (st<=chpts[j,1] && ed>=chpts[j,2]) || 
                (st>=chpts[j,1] && st<=chpts[j,2]) || 
                (ed>=chpts[j,1] && ed<=chpts[j,2]) )
            {
                overlap <- 1;
                break;
            }
        }
        if(overlap == 0){
            ret.v[i] <- TRUE
        }
            
    }

    return(ret.v)
}
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


`update.bic.baseline` <-
function(this.y,segbag,carriers,this.sum.deltasq,this.sum.deltasqwin,newchpt,samp.ord,win,y.var, 
                                PRIOR.CARRIER.PROB=0.1, PRIOR.SEGCOUNT = 10, do.plots=TRUE){
#    # for debugging:
#    this.y=t(y)          ##### remember to back transpose later.
#    newchpt=c(st,ed)
#    this.sum.deltasqwin=sum.deltasqwin[1:nrow(segbag)]
#    this.sum.deltasq = sum.deltasq[1:nrow(segbag)]
    
    
    T = nrow(this.y)
    N = ncol(this.y)
    
    m = nrow(segbag)+1    
    w = newchpt[2]-newchpt[1]+1
    
    # Each of these are length N vectors.
    if(newchpt[2]-newchpt[1]+1 >= 2){
        delta.hat = (apply(this.y[newchpt[1]:newchpt[2],],2,sum) - (w/T)*apply(this.y,2,sum))/(w*(1-w/T))
    } else {
     delta.hat = (this.y[newchpt[1]:newchpt[2],] - (w/T)*apply(this.y,2,sum))/(w*(1-w/T))
    }
    
    effect.hat = delta.hat^2/y.var
    sum.deltasq.loc = cumsum(effect.hat[samp.ord])  
    sum.deltasqwin.loc = cumsum(effect.hat[samp.ord]*w*(1-w/T))
    
    PRIOR.Nm = PRIOR.SEGCOUNT*N
    PRIOR.Carriers = PRIOR.CARRIER.PROB*PRIOR.Nm
    
    # Computing the prior probability of J, assuming a single phat=P(carrier) over all segments.
    sumJi = sum(carriers) + c(1:N) 
    Jprior.one.phat = sumJi*log((sumJi+PRIOR.Carriers)/(N*m+PRIOR.Nm)) + (N*m-sumJi)*log((N*m+PRIOR.Nm-sumJi-PRIOR.Carriers)/(N*m+PRIOR.Nm))
    
    # Computing the prior probability of J, assuming a different phat=P(carrier) over all segments.
    J = apply(carriers,2,sum)
    Jprior.sep.phat= sum( J*log((J+PRIOR.Carriers)/(N+PRIOR.Nm)) + (N-J)*log((N+PRIOR.Nm-J-PRIOR.Carriers)/(N+PRIOR.Nm)))
    Jprior.sep.phat = Jprior.sep.phat + c(1:N)*log((c(1:N)+PRIOR.Carriers)/(N+PRIOR.Nm)) + (N-c(1:N))*log((N+PRIOR.Nm-c(1:N)-PRIOR.Carriers)/(N+PRIOR.Nm))
         
    tauprior = - (lchoose(T,m) + log(m))    

    maxlik = 0.5*(sum(this.sum.deltasqwin) + sum.deltasqwin.loc)
    
    # do not take log if sum.deltasq, sum.deltasqwin is ZERO.
    
    betaprior.terms =  -sumJi/2 -(sumJi/2)*log((sum(this.sum.deltasqwin)+sum.deltasqwin.loc)/sumJi) 
    
    randomwalk.terms = 2*m*(1.6-1.5) - sum(log(this.sum.deltasq)) - log(sum.deltasq.loc)
    

    mbic.one.phat = maxlik + tauprior + Jprior.one.phat + betaprior.terms + randomwalk.terms 
    mbic.sep.phat = maxlik + tauprior + Jprior.sep.phat + betaprior.terms + randomwalk.terms 

    
    par(mfrow=c(3,1))
    plot(sqrt(effect.hat[samp.ord]*w*(1-w/T)),ylab="Chi")
    plot(maxlik,type="b", ylim=c(min(mbic.sep.phat, na.rm=TRUE),max(maxlik)))
    points(mbic.one.phat,type="b",col="red")
    segments(which.max(mbic.one.phat),0,which.max(mbic.one.phat), max(maxlik),col="red",lwd=1)
    points(mbic.sep.phat,type="b",col="blue")
    segments(which.max(mbic.sep.phat),0,which.max(mbic.sep.phat), max(maxlik),col="blue",lwd=1)
    plot(Jprior.one.phat,type="b",col="red", ylim=c(min(c(Jprior.sep.phat,betaprior.terms),na.rm=TRUE),0), ylab="Penalties",
        main=paste("PRIOR.SEGCOUNT=",PRIOR.SEGCOUNT,", PRIOR.CARRIER.PROB=",PRIOR.CARRIER.PROB))
    points(Jprior.sep.phat,type="b",col="blue")
    points(betaprior.terms,type="b",col="darkgreen")
    points(randomwalk.terms,type="b",col="orange")
     
    list(mbic.sep.phat = mbic.sep.phat, mbic.one.phat = mbic.one.phat, sum.deltasq.loc = sum.deltasq.loc, sum.deltasqwin.loc = sum.deltasqwin.loc)
}



`vu` <-
function(x,maxn=1000,do.approx=(abs(x)<0.2)){
    x=as.matrix(x)
    vux = matrix(0,nrow=nrow(x),ncol=ncol(x))

#    if(sum(is.na(do.approx)) > 0 && msscan.debug.trace){
#        #check for errors here
#        print(do.approx)
#        scan()
#        print(x)
#        scan()
#    }


    if(is.logical(do.approx)){
        if(sum(do.approx) > 0)   vux[do.approx] = exp(-x[do.approx]*0.583);
    }else if(length(do.approx) > 0){
        vux[do.approx] = exp(-x[do.approx]*0.583);
    }
  

    if (sum(do.approx)<length(x)){
        notdo.approx.ix = which(!do.approx)
        n=matrix(c(1:maxn),nrow=1,ncol=maxn)
        summands = pnorm(-0.5*matrix(x[notdo.approx.ix],nrow=length(notdo.approx.ix),ncol=1)%*%sqrt(n))/matrix(rep(n,length(notdo.approx.ix)),ncol=length(n), byrow=TRUE)
        expterm = -2*apply(summands,1,"sum");
        vux[notdo.approx.ix] = (2/x[notdo.approx.ix]^2)*exp(expterm);
    }

    vux
}

`wpca` <-
function(x, w, alpha.init =NULL,  r.init = NULL, max.iter=100, epsilon=0.01, x.lab=NULL, plots=TRUE){
    n = nrow(x)
    m = ncol(x)
    
    if(is.null(r.init)){
        r = rep(1,m)
        alpha = rep(0,m)
    } else {
        r = r.init
        alpha = alpha.init
    }

    x.vec = matrix(data=x, nrow=m*n, ncol=1)
    w.vec = matrix(data=1/w, nrow=m*n, ncol=1)
    x.hat = matrix(data=0, nrow=n, ncol=m)
    
    I = matrix(data=0, nrow=m*n, ncol=m)
    for(i in 0:(m-1)){
        I[i*n+c(1:n),i+1] = 1
    }
    
    r.prev=rep(0,m)
    iter=0
    converged=FALSE
    
    while(TRUE){
        iter=iter+1
        # estimate u given r.
        
        u.pred = matrix(ncol=n, nrow=0)
        for(i in 1:m){
            u.pred = rbind(u.pred, r[i]*diag(n))
        } 
        
        y = x.vec - I%*%alpha  # take out the platform specific baselines.
        lm1<-lm(y~u.pred-1, weights=w.vec)
        u = lm1$coefficients
    
        # estimate r, alpha given u.
        
        r.pred =  matrix(data=0, nrow=(m-1)*n, ncol=(m-1)*2)
        for(i in 0:(m-2)){
            r.pred[i*n+c(1:n),i+1] = 1
            r.pred[i*n+c(1:n),m-1+i+1] = u
        } 
        y = x.vec[1:(n*(m-1))]   # use only the first (m-1) platforms.
        lm2<-lm(y~r.pred-1, weights=w.vec[1:(n*(m-1))])
        alpha[1:(m-1)] = lm2$coefficients[1:(m-1)]
        r[1:(m-1)] = lm2$coefficients[m:(2*(m-1))]
        alpha[m] = mean(x.vec[(n*(m-1)+1):(n*m)] - u)
        r[m]=1
      
        if(sqrt(crossprod(r-r.prev)/m)< epsilon){ 
            converged=TRUE
            break
        }
        if(iter > max.iter){
            cat("Weighted PCA stopped at ",max.iter," iterations.\n")
            break
        }
        
        r.prev = r
        
        cat("r: ", r, "\n")
        
        # ----------------- diagnostic plot --------------------
        
        for(i in 1:m){
            x.hat[,i] = alpha[i] + r[i]*u
        }
        sd=sqrt(1/w)
        
        
        if(plots){
            if(m>2){ 
                par(mfrow=c(m-1,m-1))
            } else{
                par(mfrow=c(1,2))
            }
            for(i in 1:(m-1)){
                for(j in (i+1):m){
                   scatter.ci(x[,i], x[,j], 2*sd[,i], 2*sd[,j], cex=1.5, xlab=x.lab[[i]], ylab=x.lab[[j]], cex.lab=1.5)
                   abline(alpha[j]-alpha[i]*r[j]/r[i], r[j]/r[i])
                   
                   points(x.hat[,i], x.hat[,j], col="red")
                 #  for(a in 1:n){
                 #   segments(x[a,i],x[a,j], x.hat[a,i], x.hat[a,j], col="red")
                 #  }     
                }
            } 
            
            # specific to Illumina, Affymetrix, Agilent.
            col=c(1:m)
            pch=c(17:(17+m-1))
            
            plot(x[,1],type="p",col=col[1], ylim=c(min(x),max(x)+0.1), cex=2, pch=pch[1], xlab="Region index", ylab="Fitted values", cex.lab=1.5)
            lines(x.hat[,1], col=col[1])
            for(i in 2:m){
                points(x[,i],type="p",col=col[i], cex=1.5, pch=pch[i])
                lines(x.hat[,2], col=col[i])    
            }
            legend(x="topright", pch=pch, col=col, legend=x.lab, cex=1.5)
        }
    }
    
    for(i in 1:m){
        x.hat[,i] = alpha[i] + r[i]*u
    }
    
    list(fitted=x.hat, alpha=alpha, r=r, u=u, converged=converged)
}





`scatter.ci` <-
function(x,y,x.ci, y.ci, line.col="darkgray", col="black", ...){
    plot(x, y, col=col,...)
    for(i in 1:length(x)){
        segments(x[i]-x.ci[i], y[i], x[i]+x.ci[i], y[i], col=line.col)
        segments(x[i], y[i]-y.ci[i],  x[i], y[i]+y.ci[i], col=line.col)
    }
    points(x, y, col=col,...)
}






`scatter.cpfit` <-
function(x,w,x.hat,r,alpha, platform.names){
  K=ncol(x)
  par(mfrow=c(K-1,K-1))
        for(i in 1:(K-1)){
            for(j in (i+1):K){
               scatter.ci(x[,i], x[,j], 2*sd[,i], 2*sd[,j], cex=1.5, xlab=platform.names[[i]], ylab=platform.names[[j]], cex.lab=1.5)
               abline(alpha[j]-alpha[i]*r[j]/r[i], r[j]/r[i])
               
               points(x.hat[,i], x.hat[,j], col="red")
             #  for(a in 1:n){
             #   segments(x[a,i],x[a,j], x.hat[a,i], x.hat[a,j], col="red")
             #  }
            }
        }    
}










