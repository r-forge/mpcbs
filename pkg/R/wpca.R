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

