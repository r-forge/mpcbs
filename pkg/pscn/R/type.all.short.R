`type.all.short` <-
function(myname){
  typelist = c()
  nsample = dim(myname)[1]
  for (k in 1:nsample){
    a = as.character(myname[k,1])
    temp = scan(paste(a,".shortlist.txt",sep=""), what=character(0))
    if (length(temp)!=0){
      mydata1 = read.table(paste(a,".shortlist.txt",sep=""),header=TRUE)
      GGid = which(mydata1$Type=="Gain/Gain")
      GNid = which(mydata1$Type=="Gain/Normal")
      GLid = which(mydata1$Type=="Gain/Loss")
      NLid = which(mydata1$Type=="Normal/Loss")
      LLid = which(mydata1$Type=="Loss/Loss")
      SNP = mydata1$SNP.end - mydata1$SNP.begin+1
      templist = c(a, length(GGid),length(GNid), length(GLid), length(NLid), length(LLid), sum(SNP[GGid]), sum(SNP[GNid]), sum(SNP[GLid]), sum(SNP[NLid]), sum(SNP[LLid])) 
    }
    typelist = rbind(typelist,templist)
  }
  colnames(typelist) = c("sample","GG","GN","GL","NL","LL","GGSNP","GNSNP","GLSNP","NLSNP","LLSNP")
  return(typelist)                                                                                    
}

