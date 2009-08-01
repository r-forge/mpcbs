`ABtest` <-
function(cnvobj, alpha=0.05){
  
  N = dim(cnvobj$theta)[1]
  
  if (!is.null(cnvobj$chpts.new)){
    n = length(cnvobj$chpts.new)
    chptnew = c(1,cnvobj$chpts.new,(N+1))
  }else{
    n = length(cnvobj$chpts)
    chptnew = c(1,cnvobj$chpts,(N+1))
  }
  
  A = cnvobj$illumina$A
  B = cnvobj$illumina$B
  
  major.median = rep(0,(n+1))
  minor.median = rep(0,(n+1))
  total.median = rep(0,(n+1))
  
  major.result = matrix(0,2,n)
  minor.result = matrix(0,2,n)
  total.result = matrix(0,2,n)
  
  for (i in 1:n){
    if (i==1){
      id1old = which(cnvobj$newstate[(chptnew[i]):(chptnew[i+1]-1)] == 1)+(chptnew[i]-1)
      id2old = which(cnvobj$newstate[(chptnew[i]):(chptnew[i+1]-1)] == 2)+(chptnew[i]-1)
      id3old = which(cnvobj$newstate[(chptnew[i]):(chptnew[i+1]-1)] == 3)+(chptnew[i]-1)
      id4old = which(cnvobj$newstate[(chptnew[i]):(chptnew[i+1]-1)] == 4)+(chptnew[i]-1)
      majorold = c(A[id3old], B[id2old])
      minorold = c(A[id2old], B[id3old])
      totalold = c(A[id4old], B[id1old])
      major.median[i] = median(majorold)
      minor.median[i] = median(minorold)
      total.median[i] = median(totalold)
    }
    id1new = which(cnvobj$newstate[(chptnew[i+1]):(chptnew[i+2]-1)] == 1)+(chptnew[i+1]-1)
    id2new = which(cnvobj$newstate[(chptnew[i+1]):(chptnew[i+2]-1)] == 2)+(chptnew[i+1]-1)
    id3new = which(cnvobj$newstate[(chptnew[i+1]):(chptnew[i+2]-1)] == 3)+(chptnew[i+1]-1)
    id4new = which(cnvobj$newstate[(chptnew[i+1]):(chptnew[i+2]-1)] == 4)+(chptnew[i+1]-1)
    majornew = c(A[id3new], B[id2new])
    minornew = c(A[id2new], B[id3new])
    totalnew = c(A[id4new], B[id1new])
    major.median[i+1] = median(majornew)
    minor.median[i+1] = median(minornew)
    total.median[i+1] = median(totalnew)
    if(length(majorold)>0 && length(majornew)>0){major.result[,i] = WRank(majorold, majornew, alpha=alpha)}
    if(length(minorold)>0 && length(minornew)>0){minor.result[,i] = WRank(minorold, minornew, alpha=alpha)}
    if(length(totalold)>0 && length(totalnew)>0){total.result[,i] = WRank(totalold, totalnew, alpha=alpha)}
    majorold = majornew
    minorold = minornew
    totalold = totalnew
  }
  list(major.median=major.median, minor.median=minor.median, total.median=total.median, major.result=major.result, minor.result=minor.result, total.result=total.result)
}

