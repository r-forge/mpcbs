`feature.matrix.cbs` <-
function(yhat,chpts){
  T = nrow(yhat)
  N= ncol(yhat)
  if(chpts[1] != 1) chpts = c(1,chpts,T)

  feature.matrix = matrix(0,nrow=N,ncol=length(chpts)-1)
  regions = matrix(nrow=length(chpts)-1, ncol=2,data=0)
  for(i in 1:(length(chpts)-1)){
    feature.matrix[,i] = yhat[chpts[i],]
    regions[i,] = c(chpts[i], chpts[i+1])
  }
  
  list(feature.matrix=feature.matrix, feature.regions=regions)
}

