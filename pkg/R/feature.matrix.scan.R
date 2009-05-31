`feature.matrix.scan` <-
function(yhat,segments){
    pos = segments[,1]
    feature.matrix = yhat[,pos]
    feature.matrix
}

