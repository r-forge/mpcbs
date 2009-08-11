`cnvlist4` <-
function(samplename, chrs=1:23, talpha=0.01, cutvalue=1e-10, pi1=0.5, mu1.initial=1.1, mu2.initial=0.9,sigma1.initial=1, sigma2.initial=1, total0=1.5, total1 = 2.5, MIN.SNPS=20){
# total0=1.5; total1 = 2.5; level0 = 0.8; level1 = 1.2; cut1=1.9; cut2=2.1; talpha=0.01; cutvalue=1e-10; pi1=0.5; pi2=0.5; mu1.initial=1.1; mu2.initial=0.9; sigma1.initial=1; sigma2.initial=1 ; MIN.SNPS=20
# talpha: the alpha for t-test
# cutvalue: delta for iteration stop
  pi2 = 1-pi1
  if (!is.na(match(23,chrs))){
    load(paste(samplename,".Chr23.pscn.Rdata",sep=""))
    if (cnvobj$gender == "Male"){
      chrs = chrs[is.na(match(chrs,23))]
    }
  }
  # determine the normal state
  norm.total = c()
  for (chrid in chrs){
    load(paste(samplename,".Chr",chrid,".pscn.Rdata",sep=""))
    if (!is.null(cnvobj$clus)){
      n = length(cnvobj$clus)
      normid = rep(0,n) 
      
      N = length(cnvobj$newstate)
     
     chptall = c(1,cnvobj$chpts,(N+1))
      SNP.star = cnvobj$SNP.star
      if (!is.na(SNP.star)){
        if (length(which(chptall==(SNP.star+1)))==0){
          chptall = sort(c(chptall,(SNP.star+1)))
      # 3/7/2009 if add SNP.star, there may be possible that some segments have less that MIN.SNPs snps.
          throw = match(SNP.star+1,chptall)
          keepid = rep(1,length(chptall)) 
          if ((chptall[throw]-chptall[throw-1]) < MIN.SNPS){  
            keepid[throw-1] = 0
          }
          if ((chptall[throw+1]-chptall[throw]) < MIN.SNPS){
            keepid[throw+1] = 0
          }
          chptall = chptall[which(keepid==1)]            
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
        norm.total = c(norm.total,cnvobj$illumina$A[norm2],cnvobj$illumina$A[norm3], cnvobj$illumina$B[norm2],cnvobj$illumina$B[norm3])
      }
    }
  }
  mu0 = mean(norm.total)
  sigma0 = sd(norm.total)
  norm.length = length(norm.total)
  
  longlist = c()
  shortlist = c()
  
  for (chrid in chrs){
    load(paste(samplename,".Chr",chrid,".pscn.Rdata",sep=""))
    N = length(cnvobj$intensity[,1])
     
    chptall = c(1,cnvobj$chpts,(N+1))
  
    SNP.star = cnvobj$SNP.star
    if (!is.na(SNP.star)){
      if (length(which(chptall==(SNP.star+1)))==0){
        chptall = sort(c(chptall,(SNP.star+1)))
        # 3/7/2009 if add SNP.star, there may be possible that some segments have less that MIN.SNPs snps.
        throw = match(SNP.star+1,chptall)
        keepid = rep(1,length(chptall)) 
        if ((chptall[throw]-chptall[throw-1]) < MIN.SNPS){  
          keepid[throw-1] = 0
        }
        if ((chptall[throw+1]-chptall[throw]) < MIN.SNPS){
          keepid[throw+1] = 0
        }
        chptall = chptall[which(keepid==1)]     
      }
    }              
    med = ABmedian(cnvobj,chptall)
    SNP.star = cnvobj$SNP.star
    Effective.pos = cnvobj$Effective.pos            
    if (length(med$beginpt)!=length(cnvobj$clus)){
      print(paste(samplename, "Chr", chrid))
    }
    if (is.null(cnvobj$clus)){
      longtemp = c()
      shorttemp = c()
    }else{
      longtemp = c()
      n = length(chptall)-1
      count = 0
      for (k in 1:n){
        temp = cnvobj$newstate[chptall[k]:(chptall[k+1]-1)]
        id2 = which(temp==2)+chptall[k]-1
        id3 = which(temp==3)+chptall[k]-1
        n1 = length(id2)
        n2 = length(id3) 
        if (n1>0 && n2>0){
          count = count+1
          if (cnvobj$clus[count]==2 || (cnvobj$clus[count]==1 && (med$total.median[count]>total1 || med$total.median[count]<total0))){
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
     # Hao: added Mar. 29, 2009.  when there are too few Het snps, there may be only one point for one category, which will cause sigma to be 0 and delta to be NaN 
              if (!is.na(delta)){
                mu1 = mu1.new
                mu2 = mu2.new
                sigma1 = sigma1.new
                sigma2 = sigma2.new
                pi1.vec = pi1.vec.new
                pi2.vec = pi2.vec.new
              }
            }
            pvalue1 = pt((mu1-mu0)/sqrt(sigma0^2/norm.length+sum((pi1.vec[!is.na(pi1.vec)]*y[!is.na(pi1.vec)]-mu1)^2)/ny^2), (norm.length+ny-2))
            pvalue2 = pt((mu2-mu0)/sqrt(sigma0^2/norm.length+sum((pi2.vec[!is.na(pi2.vec)]*y[!is.na(pi2.vec)]-mu2)^2)/ny^2), (norm.length+ny-2))
            if (pvalue1<talpha/2){
              flag1 = -1
            }else if (pvalue1>1-talpha/2){
              flag1 = 1
            }else{
              flag1 = 0
            }
            if (pvalue2<talpha/2){
              flag2 = -1
            }else if (pvalue2>1-talpha/2){
              flag2 = 1
            }else{
              flag2 = 0
            }
                            
            if (flag1 == 1){
              if (flag2 == 1){
                Type = "Gain/Gain"
              }else if(flag2 == 0 ){
                Type = "Gain/Normal"
              }else{
                Type = "Gain/Loss"
              }
            }else if(flag1 == 0){
              if (flag2 == 1){
                Type = "Gain/Normal"
              }else if(flag2 == 0){
                Type = "NA.0.0"
              }else{
                Type = "Normal/Loss"
              }
            }else{
              if (flag2 == 1){
                Type = "Gain/Loss"
              }else if(flag2 == 0){
                Type = "Normal/Loss"
              }else{
                Type = "Loss/Loss"
              }
            }
            a = sort(c(mu1,mu2))[2]
            b = sort(c(mu1,mu2))[1]
            a = round(a,digits=3)
            b = round(b,digits=3)
            Value = paste(a,"/",b,sep="")
            SNP.begin = med$beginpt[count]
            SNP.end = med$endpt[count]
            SNP.length = SNP.end-SNP.begin+1
            SNP.percent = SNP.length/N*100
            SNP.percent = round(SNP.percent,digits=3)
            Pos.begin = cnvobj$pos[SNP.begin]
            Pos.end = cnvobj$pos[SNP.end]
            Pos.length = Pos.end-Pos.begin+1
            Pos.percent = Pos.length/Effective.pos*100
            Pos.percent = round(Pos.percent,digits=3)
            temp1 = c(samplename, chrid, Type, Value, mu0, SNP.begin, SNP.end, SNP.length, SNP.percent, Pos.begin, Pos.end, Pos.length, Pos.percent, sort(c(pvalue1,pvalue2),decreasing=TRUE))
            names(temp1)=c("Sample", "Chr", "Type","Value", "Normal.copy", "SNP.begin", "SNP.end", "SNP.length", "SNP.percent", "Pos.begin", "Pos.end", "Pos.length", "Pos.percent","major.pvalue","minor.pvalue")
            longtemp = rbind(longtemp,temp1)
          }
        }
      }
  }
    longlist = rbind(longlist,longtemp)
  }
  return(longlist)
}

