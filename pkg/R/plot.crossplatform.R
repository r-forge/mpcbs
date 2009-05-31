`plot.crossplatform` <-
function(pos,y,anchor,chpts=NULL,ranking=NULL,yhat=NULL,platform.names,num.panels=NA, xlim=NULL,...){
    K = length(y)
    
    if(is.na(num.panels)) num.panels=K
    
    par(mfrow=c(num.panels,1))
    if(is.null(xlim)){
        xlim=0
        
        xlim=c(min(anchor$merged.pos),xmax=max(anchor$merged.pos))
    } 
        
    for(k in 1:K){
        ymax = max(y[[k]]) + 0.2* max(y[[k]])
        ymin = min(y[[k]]) - 0.2*abs(min(y[[k]]))
    
        plot(pos[[k]],y[[k]], xlim=xlim,ylim=c(ymin,ymax),main=platform.names[[k]],ylab="",xlab="base position", ...)
        if(!is.null(yhat)) lines(pos[[k]],yhat[[k]], col="red",lwd=2)
        if(length(chpts)>0) {
            segments(anchor$merged.pos[chpts],rep(ymin,length(chpts)),anchor$merged.pos[chpts],rep(ymax,length(chpts)),col="black",lwd=1)
            
            if(!is.null(ranking)){
                labels = format(ranking)
                text(anchor$merged.pos[chpts], rep(max(y[[k]]),length(chpts)), labels,adj = c(0,-1))
            }
        }
    }    
}

