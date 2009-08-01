`illumina.backtransform` <-
function(cnvobj, segment=TRUE){
     logR = cnvobj$illumina$logR
     theta = cnvobj$illumina$theta
     A = cnvobj$illumina$A
     B = cnvobj$illumina$B


    if(segment && !is.null(cnvobj$chpts)){
        chpts = c(1,cnvobj$chpts,length(cnvobj)+1)
        A.fit = matrix(NA,nrow=length(chpts)-1,ncol=4)
        B.fit = matrix(NA,nrow=length(chpts)-1,ncol=4)
        logR.fit = matrix(0,nrow=length(chpts)-1,ncol=1)
        theta.fit = matrix(NA,nrow=length(chpts)-1,ncol=4)
        for(i in c(1:(length(chpts)-1))){

        
            loc = c(chpts[i]:(chpts[i+1]-1))
            A.sub = A[loc]
            B.sub = B[loc]
            theta.sub = theta[loc]
            mixstate.sub= cnvobj$mixstate[loc]
            
#            colors=c("blue","red","green","black")
#            plot(A.sub,B.sub,col=colors[mixstate.sub])
         
            logR.fit[i] = median(logR[loc])
            A.fit[i,1] = median(A.sub[mixstate.sub==1])
            B.fit[i,1] = median(B.sub[mixstate.sub==1])
            theta.fit[i,1] = median(theta.sub[mixstate.sub==1])
            A.fit[i,4] = median(A.sub[mixstate.sub==4])
            B.fit[i,4] = median(B.sub[mixstate.sub==4])
            theta.fit[i,4] = median(theta.sub[mixstate.sub==4])
            if(cnvobj$theta.seg[median(loc),2] == 0){
                A.fit[i,2] = median(A.sub[mixstate.sub==2 | mixstate.sub==3])
                B.fit[i,2] = median(B.sub[mixstate.sub==2 | mixstate.sub==3])
                theta.fit[i,2] = median(theta.sub[mixstate.sub==2 | mixstate.sub==3])
                A.fit[i,3] = A.fit[i,2]
                B.fit[i,3] = B.fit[i,2]
                theta.fit[i,3] = theta.fit[i,2]
            } else {
                A.fit[i,2] = median(A.sub[mixstate.sub==2])
                B.fit[i,2] = median(B.sub[mixstate.sub==2])
                theta.fit[i,2] = median(theta.sub[mixstate.sub==2])
                A.fit[i,3] = median(A.sub[mixstate.sub==3])
                B.fit[i,3] = median(B.sub[mixstate.sub==3])
                theta.fit[i,3] = median(theta.sub[mixstate.sub==3])

            }   
            
            
#            par(mfrow=c(2,1))
#            plot(A.sub, col=colors[mixstate.sub], pch=15)
#            for(k in c(1:4)) lines(rep(A.fit[i,k],length(loc)), col=colors[k])
#A            plot(B.sub, col=colors[mixstate.sub], pch=15)
 #           for(k in c(1:4)) lines(rep(B.fit[i,k],length(loc)), col=colors[k])
   
        }

    } 
    
    
    cnvobj$illumina$theta.fit=theta.fit
    cnvobj$illumina$logR.fit=logR.fit
    cnvobj$illumina$A.fit=A.fit
    cnvobj$illumina$B.fit=B.fit
    
    cnvobj
}

