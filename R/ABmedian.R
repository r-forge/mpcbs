`ABmedian` <-
function(cnvobj, chptall){
  n = length(chptall)-1
 
  A = cnvobj$illumina$A
  B = cnvobj$illumina$B

  major.median = c()
  minor.median = c()
  total.median = c()
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
      major.median = c(major.median,median(major))
      minor.median = c(minor.median,median(minor))
    # 3/7/2009 add the possibility that there are no homogenous in the segment
      if (length(total)>0){
        total.median = c(total.median,median(total))
      }else{
        total.meidan = c(total.median,median(major)+median(minor))
      }
      beginpt = c(beginpt,chptall[i])
      endpt = c(endpt,(chptall[i+1]-1))
    }
  }  
  list(beginpt=beginpt, endpt=endpt, major.median=major.median, minor.median=minor.median, total.median=total.median)
}

