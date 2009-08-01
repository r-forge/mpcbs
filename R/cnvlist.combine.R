`cnvlist.combine` <-
function(mydata, chrid, samplename, chrs=1:23, cutvalue=1e-10, pi1=0.5, pi2=0.5, mu1.initial=1.1, mu2.initial=0.9,sigma1.initial=1, sigma2.initial=1, total0=1.5, total1 = 2.5, MIN.SNPS=20){
# cutvalue=1e-10; pi1=0.5; pi2=0.5; mu1.initial=1.1; mu2.initial=0.9; sigma1.initial=1; sigma2.initial=1; total0=1.5; total1 = 2.5
  load(paste(samplename,".Chr",chrid,".pscn.Rdata",sep=""))
  SNP.star = cnvobj$SNP.star
  SNP.begin = mydata$SNP.begin
  SNP.end = mydata$SNP.end
  n = dim(mydata)[1]
  for (i in 1:(n-1)){
    if ( (SNP.end[i]+1) == SNP.begin[i+1] && mydata$Type[i]==mydata$Type[i+1]){
      if (is.na(SNP.star)){
        SNP.end[i] = NA
        SNP.begin[i+1] = NA
      }else if (SNP.end[i] != SNP.star){
        SNP.end[i] = NA
        SNP.begin[i+1] = NA
      }
    }
  }
  keepid = rep(1,n)
  for (i in 2:n){
    if (is.na(SNP.begin[i])){
      keepid[i] = 0
    }
  }
  recal = rep(0,n)
  for (i in 1:n){
    if (is.na(SNP.end[i])){
      recal[i] = 1
    }
  }
  recal = recal[which(keepid==1)]
  value2 = mydata$Value[which(keepid==1)]
  pvalue1.vec = mydata$major.pvalue[which(keepid==1)]
  pvalue2.vec = mydata$minor.pvalue[which(keepid==1)]
  Type.new = mydata$Type[which(keepid==1)]
  SNP.begin.new = SNP.begin[which(!is.na(SNP.begin))]
  SNP.end.new = SNP.end[which(!is.na(SNP.end))]
  Pos.begin.new = cnvobj$pos[SNP.begin.new]
  Pos.end.new = cnvobj$pos[SNP.end.new]
  m = length(SNP.begin.new)
  value.new = c()
  pvalue1.new = c()
  pvalue2.new = c()
  if (sum(recal)>0){
    if (!is.na(match(23,chrs))){
      load(paste(samplename,".Chr23.pscn.Rdata",sep=""))
      if (cnvobj$gender == "Male"){
        chrs = chrs[is.na(match(chrs,23))]
      }
    }
      # determine the normal state
    norm.total = c()
    for (ii in chrs){
      load(paste(samplename,".Chr",ii,".pscn.Rdata",sep=""))
      if (!is.null(cnvobj$clus)){
        n = length(cnvobj$clus)
        normid = rep(0,n) 
        
        N = length(cnvobj$newstate)
       
       chptall = c(1,cnvobj$chpts,(N+1))
        SNP.star = cnvobj$SNP.star
        if (!is.na(SNP.star)){
          if (length(which(chptall==(SNP.star+1)))==0){
            chptall = sort(c(chptall,(SNP.star+1)))
            throw2 = match(SNP.star + 1, chptall)
            keepid2 = rep(1, length(chptall))
            if ((chptall[throw2] - chptall[throw2 - 1]) < MIN.SNPS) {
              keepid2[throw2 - 1] = 0
            }
            if ((chptall[throw2 + 1] - chptall[throw2]) < MIN.SNPS) {
              keepid2[throw2 + 1] = 0
            }
            chptall = chptall[which(keepid2 == 1)]                         
          }  
        } 
                             
        med = ABmedian(cnvobj,chptall)
             
        for (i in 1:n){
          if (cnvobj$clus[i]==1 && med$total.median>total0 && med$total.median<total1){
            normid[i] = 1
          }
        }
        if (sum(normid)>0){
          norm2 = c()
          norm3 = c()
          for (i in 1:n){
            if (normid[i]==1){
              temp2 = which(cnvobj$newstate[med$beginpt[i]:med$endpt[i]]==2)+med$beginpt[i]-1
              temp3 = which(cnvobj$newstate[med$beginpt[i]:med$endpt[i]]==3)+med$beginpt[i]-1
              norm2 = c(norm2,temp2)
              norm3 = c(norm3,temp3)
            }
          }
          norm.total = c(norm.total,cnvobj$intensity[norm2,2],cnvobj$intensity[norm3,2], cnvobj$intensity[norm2,1],cnvobj$intensity[norm3,1])
        }
      }
    }
    mu0 = mean(norm.total)
    sigma0 = sd(norm.total)
    norm.length = length(norm.total)
  }  
  
  load(paste(samplename,".Chr",chrid,".pscn.Rdata",sep=""))
  value2 = mydata$Value[which(keepid==1)]
  pvalue1.vec = mydata$major.pvalue[which(keepid==1)]
  pvalue2.vec = mydata$minor.pvalue[which(keepid==1)]
  Type.new = mydata$Type[which(keepid==1)]
  SNP.begin.new = SNP.begin[which(!is.na(SNP.begin))]
  SNP.end.new = SNP.end[which(!is.na(SNP.end))]
  Pos.begin.new = cnvobj$pos[SNP.begin.new]
  Pos.end.new = cnvobj$pos[SNP.end.new]
  m = length(SNP.begin.new)
  value.new = c()
  pvalue1.new = c()
  pvalue2.new = c()

  for (k in 1:m){
    if (recal[k]==1){
      temp = cnvobj$newstate[SNP.begin.new[k]:SNP.end.new[k]]
      id2 = which(temp==2)+SNP.begin.new[k]-1
      id3 = which(temp==3)+SNP.begin.new[k]-1
      n1 = length(id2)
      n2 = length(id3) 
      y = c(cnvobj$illumina$A[id2],cnvobj$illumina$A[id3],cnvobj$illumina$B[id2],cnvobj$illumina$B[id3])
      ny = length(y)
      mu1 = mu1.initial
      mu2 = mu2.initial
      sigma1 = sigma1.initial
      sigma2 = sigma2.initial
      pi1.vec.new = rep(0,ny)
      pi2.vec.new = rep(0,ny)
      delta = 1
      while (!is.na(delta) && delta>cutvalue){
        for (i in 1:ny){
          pi1.vec.new[i] = exp(-(mu1-y[i])^2/(2*sigma1^2))*pi1/sigma1/(exp(-(mu1-y[i])^2/(2*sigma1^2))*pi1/sigma1+exp(-(mu2-y[i])^2/(2*sigma2^2))*pi2/sigma2)
          pi2.vec.new[i] = exp(-(mu2-y[i])^2/(2*sigma2^2))*pi2/sigma2/(exp(-(mu1-y[i])^2/(2*sigma1^2))*pi1/sigma1+exp(-(mu2-y[i])^2/(2*sigma2^2))*pi2/sigma2)
        }
        mu1.new = sum(pi1.vec.new[!is.na(pi1.vec.new)]*y[!is.na(pi1.vec.new)])/sum(pi1.vec.new[!is.na(pi1.vec.new)])
        sigma1.new = sqrt(sum(pi1.vec.new[!is.na(pi1.vec.new)]*(y[!is.na(pi1.vec.new)]-mu1.new)^2)/sum(pi1.vec.new[!is.na(pi1.vec.new)]))
        mu2.new = sum(pi2.vec.new[!is.na(pi2.vec.new)]*y[!is.na(pi2.vec.new)])/sum(pi2.vec.new[!is.na(pi2.vec.new)])
        sigma2.new = sqrt(sum(pi2.vec.new[!is.na(pi2.vec.new)]*(y[!is.na(pi2.vec.new)]-mu2.new)^2)/sum(pi2.vec.new[!is.na(pi2.vec.new)]))
        delta = (mu1.new-mu1)^2+(mu2.new-mu2)^2+(sigma1.new-sigma1)^2+(sigma2.new-sigma2)^2
        if (!is.na(delta)){
          mu1 = mu1.new
          mu2 = mu2.new
          sigma1 = sigma1.new
          sigma2 = sigma2.new
          pi1.vec = pi1.vec.new
          pi2.vec = pi2.vec.new
        }
      }
      a = sort(c(mu1,mu2))[2]
      b = sort(c(mu1,mu2))[1]
      a = round(a,digits=3)
      b = round(b,digits=3)
      Value = paste(a,"/",b,sep="") 
      value.new = c(value.new, Value)
      pvalue1 = pt((mu1-mu0)/sqrt(sigma0^2/norm.length+sum((pi1.vec[!is.na(pi1.vec)]*y[!is.na(pi1.vec)]-mu1)^2)/ny^2), (norm.length+ny-2))
      pvalue2 = pt((mu2-mu0)/sqrt(sigma0^2/norm.length+sum((pi2.vec[!is.na(pi2.vec)]*y[!is.na(pi2.vec)]-mu2)^2)/ny^2), (norm.length+ny-2))
      pvalue1.new = c(pvalue1.new, as.character(pvalue1))
      pvalue2.new = c(pvalue2.new, as.character(pvalue2))
    }else{
      value.new = c(value.new, as.character(value2[k]))
      pvalue1.new = c(pvalue1.new, (pvalue1.vec[k]))
      pvalue2.new = c(pvalue2.new, (pvalue2.vec[k]))
    }
  }
  mydata.new = cbind(rep(samplename,m), rep(chrid,m), as.character(Type.new), value.new, rep(mydata$Normal.copy[1],m), SNP.begin.new, SNP.end.new, Pos.begin.new, Pos.end.new,pvalue1.new, pvalue2.new)
  colnames(mydata.new) = c("Sample", "Chr", "Type", "Value", "Normal.copy", "SNP.begin", "SNP.end", "Pos.begin", "Pos.end", "major.pvalue", "minor.pvalue")
  return(mydata.new)
}

