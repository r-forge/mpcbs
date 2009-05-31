`matrix.max` <-
function(Z){
    # If Z has NA values treat as -Inf
    na.inds = which(is.na(Z),arr.ind=TRUE)
    Z[na.inds]=-Inf

    row.max.col=max.col(Z, ties.method=c("random", "first", "last"))
    max.row = which.max(Z[cbind(1:nrow(Z),row.max.col)])
    c(max.row,row.max.col[max.row])
}

