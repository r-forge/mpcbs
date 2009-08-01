`plot.cnv2` <-
function(mydata, mu0, n, Ntest, samplename, chrid, alpha=0.01, gaincol="red", normalcol="green", losscol="blue", elsecol = "grey", linewidth = 3, threshold = 0.1){
# alpha=0.01; gaincol="red"; normalcol="green"; losscol="blue"; elsecol = "grey"; linewidth = 3
    load(paste(samplename,".Chr",chrid,".pscn.Rdata",sep=""))
    N = length(cnvobj$pos)
    x = 1:N
    major.pvalue = mydata$major.pvalue
    minor.pvalue = mydata$minor.pvalue                                                                                                                 
    major = c()
    minor = c()
    for (j in 1:n){
      if (major.pvalue[j]<(1-alpha/(Ntest*2)) && major.pvalue[j]>(alpha/(Ntest*2))){
        if (minor.pvalue[j]<(1-alpha/(Ntest*2)) && minor.pvalue[j]>(alpha/(Ntest*2))){
          valueM = as.numeric(strsplit(as.character(mydata$Value[j]),"/")[[1]][1])
          major = c(major,valueM)
          valuem = as.numeric(strsplit(as.character(mydata$Value[j]),"/")[[1]][2])
          minor = c(minor,valuem)
        }else{
          major = c(major,mu0)
          valuem = as.numeric(strsplit(as.character(mydata$Value[j]),"/")[[1]][2])
          minor = c(minor,valuem)
        }
      }else if(major.pvalue[j]>(1-alpha/(Ntest*2))){
        valueM = as.numeric(strsplit(as.character(mydata$Value[j]),"/")[[1]][1])
        major = c(major,valueM)
        if(minor.pvalue[j]<(1-alpha/(Ntest*2)) && minor.pvalue[j]>(alpha/(Ntest*2))){
          minor = c(minor,mu0)
        }else{
          valuem = as.numeric(strsplit(as.character(mydata$Value[j]),"/")[[1]][2])
          minor = c(minor,valuem)
        }
      }else{
        valueM = as.numeric(strsplit(as.character(mydata$Value[j]),"/")[[1]][1])
        major = c(major,valueM)
        valuem = as.numeric(strsplit(as.character(mydata$Value[j]),"/")[[1]][2])
        minor = c(minor,valuem)
      }
    }
    normal = rep(mu0,N)
    for (j in 1:n){
      if (major[j]!=mu0 && minor[j]!=mu0){
        SNPid = mydata$SNP.begin[j]:mydata$SNP.end[j]
        normal[SNPid] = NA
      }
    }

    plot(x,normal,type="l",col=normalcol,lwd=linewidth,xlab="SNP index",ylab="Copy Number",main=paste(samplename," Chr",chrid),ylim=c(0,max(2,max(major))))
    for (j in 1:n){
      SNPid = mydata$SNP.begin[j]:mydata$SNP.end[j]
      if (major.pvalue[j]>(1-alpha/(Ntest*2))){
        points(SNPid,rep(major[j],length(SNPid)),type="l",col=gaincol,lwd=linewidth)
        points(rep(mydata$SNP.begin[j],2),c(mu0,major[j]),type="l",col=gaincol)
        points(rep(mydata$SNP.end[j],2),c(mu0,major[j]),type="l",col=gaincol)
      }else if (major.pvalue[j]<alpha/(Ntest*2)){
        points(SNPid,rep(major[j],length(SNPid)),type="l",col=losscol,lwd=linewidth)
        points(rep(mydata$SNP.begin[j],2),c(major[j],mu0),type="l",col=losscol)
        points(rep(mydata$SNP.end[j],2),c(major[j],mu0),type="l",col=losscol)
      }
      if (minor.pvalue[j]>(1-alpha/(Ntest*2))){
        points(SNPid,rep(minor[j],length(SNPid)),type="l",col=gaincol,lwd=linewidth)
        points(rep(mydata$SNP.begin[j],2),c(mu0,minor[j]),type="l",col=gaincol)
        points(rep(mydata$SNP.end[j],2),c(mu0,minor[j]),type="l",col=gaincol)
      }else if (minor.pvalue[j]<alpha/(Ntest*2)){ 
        points(SNPid,rep(minor[j],length(SNPid)),type="l",col=losscol,lwd=linewidth)
        points(rep(mydata$SNP.begin[j],2),c(minor[j],mu0),type="l",col=losscol)
        points(rep(mydata$SNP.end[j],2),c(minor[j],mu0),type="l",col=losscol)
      }
      if (major.pvalue[j]<(1-alpha/(Ntest*2)) && major.pvalue[j]>(alpha/(Ntest*2))&& minor.pvalue[j]<(1-alpha/(Ntest*2)) && minor.pvalue[j]>(alpha/(Ntest*2))){
        points(SNPid,rep(major[j],length(SNPid)),type="l",col=elsecol,lwd=linewidth)
        points(rep(mydata$SNP.begin[j],2),c(mu0,major[j]),type="l",col=elsecol)
        points(rep(mydata$SNP.end[j],2),c(mu0,major[j]),type="l",col=elsecol)
        points(SNPid,rep(minor[j],length(SNPid)),type="l",col=elsecol,lwd=linewidth)
        points(rep(mydata$SNP.begin[j],2),c(minor[j],mu0),type="l",col=elsecol)
        points(rep(mydata$SNP.end[j],2),c(minor[j],mu0),type="l",col=elsecol)        
      }
      if (major.pvalue[j]>(1-alpha/(Ntest*2)) && minor.pvalue[j]>(1-alpha/(Ntest*2))){
        if ((minor[j]-mu0)<threshold){
          points(SNPid,rep(minor[j],length(SNPid)),type="l",col=normalcol,lwd=linewidth)
        }
      }
      if (major.pvalue[j]<alpha/(Ntest*2) && minor.pvalue[j]<alpha/(Ntest*2)){
        if ((mu0-major[j])<threshold){
          points(SNPid,rep(major[j],length(SNPid)),type="l",col=normalcol,lwd=linewidth)
        }
      }
      if (major.pvalue[j]>(1-alpha/(Ntest*2)) && minor.pvalue[j]<alpha/(Ntest*2)){
        a = major[j]-mu0
        b = mu0-minor[j]
        if (a>b && b<threshold){
          points(SNPid,rep(minor[j],length(SNPid)),type="l",col=normalcol,lwd=linewidth)
          points(rep(mydata$SNP.begin[j],2),c(minor[j],mu0),type="l",col=normalcol)
          points(rep(mydata$SNP.end[j],2),c(minor[j],mu0),type="l",col=normalcol)
        }else if (b>a && a<threshold){
          points(SNPid,rep(major[j],length(SNPid)),type="l",col=normalcol,lwd=linewidth)
          points(rep(mydata$SNP.begin[j],2),c(mu0,major[j]),type="l",col=normalcol)
          points(rep(mydata$SNP.end[j],2),c(mu0,major[j]),type="l",col=normalcol)
        }
      }
  }
}

