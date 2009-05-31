`scatter.ci` <-
function(x,y,x.ci, y.ci, line.col="darkgray", col="black", ...){
    plot(x, y, col=col,...)
    for(i in 1:length(x)){
        segments(x[i]-x.ci[i], y[i], x[i]+x.ci[i], y[i], col=line.col)
        segments(x[i], y[i]-y.ci[i],  x[i], y[i]+y.ci[i], col=line.col)
    }
    points(x, y, col=col,...)
}

