`regroup.theta` <-
function(cnvobj,level0=0.1,level1=0.9, percent=0.15){

  if (is.null(cnvobj$sig.seg)){
        cnvobj<-segment(cnvobj)
  }
  
  N = dim(cnvobj$theta)[1]
  
  if (!is.null(cnvobj$chpts.new)){
    if (length(cnvobj$chpts.new)>0){
      n = length(cnvobj$chpts.new)+1
      chptall = c(1,cnvobj$chpts.new,(N+1))
    }else{
      n = 1
      chptall = c(1,(N+1))
    }
  }else{
    if (is.null(cnvobj$chpts)){
      n = 1
      chptall = c(0,N)
    }else{
      n = length(cnvobj$chpts)+1
      chptall = c(1,cnvobj$chpts,(N+1))
    }
  }
  
  theta = cnvobj$illumina$theta
  newstate = c()

  for (i in 1:n){
    temp = theta[chptall[i]:(chptall[i+1]-1)]
    temp2 = kmeans(temp,centers=c(1,0.9,0.1,0),algorithm="Lloyd")
    old = temp2$cluster
    if (temp2$size[2]<ceiling(length(temp)*percent/2)){
      id1 = which(temp2$cluster==1)
      temp3 = temp[id1]
      if (length(temp3)>3 && (max(temp3) - min(temp3)) > 0) {
        temp4 = kmeans(temp3,centers=2)
        if (temp4$centers[1]>temp4$centers[2]){
          old[id1] = temp4$cluster
        }else{
          old[id1] = 3-temp4$cluster
        }
      }
    }
    if (temp2$size[3]<ceiling(length(temp)*percent/2)){
      id4 = which(temp2$cluster==4)
      temp5 = temp[id4]
      if (length(temp5)>3 && (max(temp5) - min(temp5)) > 0) {
        temp6 = kmeans(temp5,centers=2)
        if (temp6$centers[1]>temp6$centers[2]){
          old[id4] = temp6$cluster+2
        }else{
          old[id4] = 5-temp6$cluster
        }
      }
    }
    newstate = c(newstate, old)
  }    
  cnvobj$newstate = newstate
  cnvobj
}

