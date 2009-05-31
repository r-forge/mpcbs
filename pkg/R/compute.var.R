`compute.var` <-
function(y, use.mean=TRUE){

    if(is.null(dim(y))){
         y = as.matrix(y,nrow=length(y),ncol=1)
    }

    N=dim(y)[2]
    T=dim(y)[1]
    y.diff<-y[1:(T-1),] - y[2:T,]
    
    if(!use.mean){
        if(N>1){
            y.var <- apply(y.diff^2,2,median)
        } else {
            y.var = median(y.diff^2)
        }
    } else {
        if(N>1){
            y.var <- apply(y.diff^2,2,mean)/2
        } else {
            y.var = mean(y.diff^2)/2
        }
    }
    
    y.var
}

