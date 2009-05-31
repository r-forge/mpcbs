`Remove.Overlap.C` <-
function(chpts){
    res1 <- .Call("RemoveOverlap", chpts[,1], chpts[,2], PACKAGE="mpcbs");
    
    i2 <- rep(FALSE, nrow(chpts))
    i2[res1 == 1] <- TRUE
    return(i2)
}

