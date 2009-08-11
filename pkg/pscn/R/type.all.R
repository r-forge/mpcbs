`type.all` <-
function(myname){
  typelist = c()
  nsample = dim(myname)[1]
  for (k in 1:nsample){
    a = as.character(myname[k,1])
    temp = scan(paste(a,".longlist.update.txt",sep=""), what=character(0))
    if (length(temp)!=0){
      mydata1 = read.table(paste(a,".longlist.update.txt",sep=""),header=TRUE)
      templist = rep(0,10)
      for (i in 1:23){
        id1 = which(mydata1$Chr==i)
        n1 = length(id1)
        if (n1>0){
          cnvlist1 = mydata1[id1,]
          temp = type.sub(cnvlist1,cnvlist1$Normal.copy[1],n1,dim(mydata1)[1],a,i)
          templist = templist+temp 
        }
      }
      typelist = rbind(typelist,c(a,templist))
    }
  }
  colnames(typelist) = c("sample","GG","GN","GL","NL","LL","GGSNP","GNSNP","GLSNP","NLSNP","LLSNP")
  return(typelist)
}

