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

