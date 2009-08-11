`group.number` <-
function(samplename, chrs = 1:23, gender=NULL, MIN.SNPS = 20){
# dataflag is used to treat the original 63 samples different from the later samples because the original 63 samples use a different system in the name. dataflag=0 for the original 63 samples, while 1 for the later samples.

# determine gender first
    if (!is.na(match(23,chrs))){
      load(paste(samplename,".Chr",23,".Segment.Rdata",sep=""))
      id2 = length(which(cnvobj$newstate==2))
      id3 = length(which(cnvobj$newstate==3))
      N = length(cnvobj$newstate)
      percent = (id2+id3)/N
      if (percent>0.1){
        gender = "Female"
      }else{
        gender = "Male"
      }
    }
    
   for (jj in chrs){
    load(paste(samplename,".Chr",jj,".Segment.Rdata",sep=""))
    N = length(cnvobj$pos)
    if (is.null(cnvobj$chpts)){
      chptall = c(1,(N+1))
    }else{
      if (length(cnvobj$chpts.new)==0){
        chptall = c(1,(N+1))
      }else{
        chptall = c(1,cnvobj$chpts.new,(N+1))
      }
    }
    temp = centro(cnvobj)
    SNP.star = temp[4]
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
    clus = clustertest(cnvobj,chptall)
    cnvobj$SNP.star = SNP.star
    cnvobj$Effective.pos = temp[3]
    cnvobj$gender = gender
    cnvobj$clus = clus[5,]
    cnvobj$chpts = cnvobj$chpts.new
    cnvobj$chpts.new = NULL
    save(cnvobj,file=paste(cnvobj$label,".pscn.Rdata",sep=""))
  }
}

