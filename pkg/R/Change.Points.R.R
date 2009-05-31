`Change.Points.R` <-
function(ix, T, lookup){

    from.row <- (ix[1:lookup] %% T)
    from.col <- as.integer( ix[1:lookup] / T) + 1
    
    v1 <- from.row + 1
    v2 <- from.col + from.row
    return(cbind(v1, v2))
}

