`ABmean` <-
function(cnvobj, chptall){
  n = length(chptall)-1
 
  A = cnvobj$illumina$A
  B = cnvobj$illumina$B

  major.mean = c()
  minor.mean = c()
  total.mean = c()
  beginpt = c()
  endpt = c() 
  
  for (i in 1:n){
    temp = cnvobj$newstate[chptall[i]:(chptall[i+1]-1)]
    id2 = which(temp==2)+chptall[i]-1
    id3 = which(temp==3)+chptall[i]-1
    n1 = length(id2)
    n2 = length(id3)
    if (n1>0 && n2>0){
      major = c(A[id3],B[id2])
      minor = c(A[id2],B[id3])
      id4 = which(temp==4)+chptall[i]-1
      id1 = which(temp==1)+chptall[i]-1
      total = c(A[id4],B[id1])
      major.mean = c(major.mean,mean(major))
      minor.mean = c(minor.mean,mean(minor))
      total.mean = c(total.mean,mean(total))
      beginpt = c(beginpt,chptall[i])
      endpt = c(endpt,(chptall[i+1]-1))
    }
  }  
  list(beginpt=beginpt, endpt=endpt, major.mean=major.mean, minor.mean=minor.mean, total.mean=total.mean)
}

