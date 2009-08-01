`pscnlist` <-
function(samplename, chrs=1:23, MIN.SNPS=20, alpha=0.01, FWalpha=0.01, threshold=0.18, hard.threshold=0.125, ...){
  longlist = cnvlist4(samplename=samplename, chrs=chrs, MIN.SNPS=20, talpha=alpha, ...)
  write.table(longlist, file=paste(samplename,".longlist.txt",sep=""),row.names=FALSE, sep="\t", quote=FALSE)
  mydata = read.table(paste(samplename,".longlist.txt",sep=""),header=TRUE)
  longlist2 = cnvlist.update(mydata, alpha=FWalpha, threshold=threshold, hard.threshold=hard.threshold)
  longid = which(!is.na(longlist2$Type))
  longlist2 = longlist2[longid,]
  write.table(longlist2, file=paste(samplename,".longlist.update.txt",sep=""),row.names=FALSE, sep="\t", quote=FALSE)
  mydata2 = read.table(paste(samplename,".longlist.update.txt",sep=""),header=TRUE)
  mydata2 = mydata2[which(!is.na(mydata2$Type)),]
  shortlist = c()
  for (i in chrs){
    id = which(mydata2$Chr==i)
    if (length(id)>0){
      chrdata = mydata2[id,]
      if (length(id)==1){
        temp = chrdata[c(1:7,10:11, 14:15)]
        temp$Chr = as.character(temp$Chr)
        temp$SNP.begin = as.character(temp$SNP.begin)
        temp$SNP.end = as.character(temp$SNP.end)
        temp$Pos.begin = as.character(temp$Pos.begin)
        temp$Pos.end = as.character(temp$Pos.end)
        temp$major.pvalue = as.character(temp$major.pvalue)
        temp$minor.pvalue = as.character(temp$minor.pvalue)
        shortlist = rbind(shortlist, temp )
      }else{
        chrdata2 = cnvlist.combine(chrdata, chrid=i, samplename=samplename, chrs=chrs, MIN.SNPS=MIN.SNPS, ...)
        shortlist = rbind(shortlist,chrdata2)
      }
    }
  }
  write.table(shortlist, file=paste(samplename,".shortlist.txt",sep=""),row.names=FALSE, sep="\t", quote=FALSE)      
}

