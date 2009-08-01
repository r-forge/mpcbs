`type.sub` <-
function(mydata, mu0, n, Ntest, samplename, chrid, alpha=0.01, threshold = 0.15){
  load(paste(samplename,".Chr",chrid,".pscn.Rdata",sep=""))
  major.pvalue = mydata$major.pvalue
  minor.pvalue = mydata$minor.pvalue
  GG = 0
  GN = 0
  GL = 0
  NL = 0
  LL = 0
  GGSNP = 0
  GNSNP = 0
  GLSNP = 0
  NLSNP = 0
  LLSNP = 0
  for (j in 1:n){
    if (major.pvalue[j]>(1-alpha/(Ntest*2))){
      if (minor.pvalue[j]>(1-alpha/(Ntest*2)) && (as.numeric(strsplit(as.character(mydata$Value[j]),"/")[[1]][2])-mu0)>threshold){
        GG = GG+1
        GGSNP = GGSNP + mydata$SNP.length[j]
      }else if (minor.pvalue[j]<alpha/(Ntest*2)){
        a = as.numeric(strsplit(as.character(mydata$Value[j]),"/")[[1]][1])-mu0
        b = mu0-as.numeric(strsplit(as.character(mydata$Value[j]),"/")[[1]][2])
        if (a>b && b<threshold){
          GN = GN+1
          GNSNP = GNSNP + mydata$SNP.length[j]
        }else if (b>a && a<threshold){
          NL = NL+1
          NLSNP = NLSNP + mydata$SNP.length[j]
        }else{
          GL = GL+1
          GLSNP = GLSNP + mydata$SNP.length[j]
        }
      }else{                                                                                                                                                                            
        GN = GN+1
        GNSNP = GNSNP + mydata$SNP.length[j]
      }
    }else if (major.pvalue[j]<alpha/(Ntest*2)){
      if (minor.pvalue[j]<alpha/(Ntest*2) && (mu0-as.numeric(strsplit(as.character(mydata$Value[j]),"/")[[1]][1]))>threshold ){
        LL = LL+1
        LLSNP = LLSNP + mydata$SNP.length[j]
      }else if (minor.pvalue[j]<alpha/(Ntest*2) && (mu0-as.numeric(strsplit(as.character(mydata$Value[j]),"/")[[1]][1]))<threshold ){
        NL = NL+1
        NLSNP = NLSNP + mydata$SNP.length[j]
      }
    }else{
      if (minor.pvalue[j]<alpha/(Ntest*2)){
        NL = NL+1
        NLSNP = NLSNP + mydata$SNP.length[j]
      }
    }
  }
  return(c(GG,GN,GL,NL,LL,GGSNP,GNSNP,GLSNP,NLSNP,LLSNP))
}

