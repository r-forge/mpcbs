`clustertest` <-
function(cnvobj,chptall){
  n = length(chptall)-1
    
  test1 = c()
  test2 = c()
  beginpt = c()
  endpt = c()
  result = c()

  for (i in 1:n){
    temp = cnvobj$newstate[chptall[i]:(chptall[i+1]-1)]
    id2 = which(temp==2)+chptall[i]-1
    id3 = which(temp==3)+chptall[i]-1
    mat1 = cnvobj$intensity[id2,]
    mat2 = cnvobj$intensity[id3,]
    n1 = length(id2)
    n2 = length(id3)
    if (n1>2 && n2>2){
      beginpt = c(beginpt,chptall[i])
      endpt = c(endpt,(chptall[i+1]-1))     
      test1 = c(test1,cluster.1(mat1,mat2))
      test2 = c(test2,cluster.2(mat1,mat2))
      m = length(test1)
      # add "as.double" on Apr. 4, 2009, otherwise may cause problem.
      if (as.double(test1[m])>as.double(test2[m])){
        result = c(result,1)
      }else{
        result = c(result,2)
      }
    }else if (n1>0 && n2>0){
  # 3/7 to be consistent with ABmedian
      beginpt = c(beginpt,chptall[i])
      endpt = c(endpt,(chptall[i+1]-1))
      test1 = c(test1,"NA")
      test2 = c(test2,"NA")
      result = c(result,1)
    }         
  }    
  return(rbind(beginpt,endpt,test1,test2,result))
}

