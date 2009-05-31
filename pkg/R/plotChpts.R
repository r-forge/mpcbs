`plotChpts` <-
function(msscan.retval){

    plot(   c(0, ncol(msscan.retval$yhat)),
        c(0, nrow(msscan.retval$yhat)),
        xlab="sample index",
        ylab="SNP index",
        col="white",
        main="")
    
    for(j in 1:nrow(msscan.retval$yhat)){
        v1 <- msscan.retval$yhat[j,-1]
        v2 <- msscan.retval$yhat[j,-ncol(msscan.retval$yhat)]
        
        i2 <- which(v2 != v1)

        if(length(i2) == 0) next

        for(k in i2){
            lines(c(k, k+1, k+1, k,k), c(j,j,j+1,j+1,j), lty="solid")
        }
    }
        
}

