`plot.cnv3` <-
function(mydata, chrid, pos=FALSE, gaincol="red", mainlab=NULL, normalcol="green", losscol="blue", elsecol = "grey", linewidth = 3, locs=NULL){
    mydata = mydata[which(mydata$Chr==chrid),]
# alpha=0.01; gaincol="red"; normalcol="green"; losscol="blue"; elsecol = "grey"; linewidth = 3
    load(paste(mydata$Sample[1],".Chr",chrid,".pscn.Rdata",sep=""))
    mu0 = mydata$Normal.copy[1]
    n = dim(mydata)[1]
    N = length(cnvobj$pos)
    Type = mydata$Type
    if (is.null(locs)){
      a = 1
      b = N
    }else{
      locs = sort(locs)
      a = max(locs[1],1)
      b = min(locs[length(locs)],N)
    }
    x = a:b
    id = rep(0,n)
    for (j in 1:n){
      if (mydata$SNP.end[j]>=a && mydata$SNP.begin[j]<=b){
        id[j] = 1
      }
    }
    mydata = mydata[which(id==1),]
    n = dim(mydata)[1]       
    major = c()
    minor = c()
    if (n > 0){
      for (j in 1:n){
        if (Type[j]=="Gain/Loss" || Type[j]=="Gain/Gain" || Type[j]=="Loss/Loss"){
          valueM = as.numeric(strsplit(as.character(mydata$Value[j]),"/")[[1]][1])
          major = c(major,valueM)
          valuem = as.numeric(strsplit(as.character(mydata$Value[j]),"/")[[1]][2])
          minor = c(minor,valuem)
        }else if (Type[j] == "Normal/Loss"){
          major = c(major,mu0)
          valuem = as.numeric(strsplit(as.character(mydata$Value[j]),"/")[[1]][2])
          minor = c(minor,valuem)
        }else if(Type[j]=="Gain/Normal"){
          valueM = as.numeric(strsplit(as.character(mydata$Value[j]),"/")[[1]][1])
          major = c(major,valueM)
          minor = c(minor,mu0)
        }  
      }
    }
    normal = rep(mu0,(b-a+1))
    if (n >0){
      for (j in 1:n){
        if (major[j]!=mu0 && minor[j]!=mu0){
          SNPid = max(1,(mydata$SNP.begin[j]-a+1)):min((mydata$SNP.end[j]-a+1),(b-a+1))
          normal[SNPid] = NA
        }
      }
    }
    if (pos){
      if (!is.null(cnvobj$SNP.star) && cnvobj$SNP.star>=a && cnvobj$SNP.star<=b){
        normal[cnvobj$SNP.star-a+1]=NA
      }
      plot(cnvobj$pos[x],normal,type="l",col=normalcol,lwd=linewidth,xlab="Position (bp)",ylab="Copy Number",main=mainlab,ylim=c(0,max(2,max(major))))
      if (n>0){
        for (j in 1:n){                                                            
          SNPid = max(a,mydata$SNP.begin[j]):min(b,mydata$SNP.end[j])
          if (Type[j] == "Gain/Gain" || Type[j]=="Gain/Normal" || Type[j]=="Gain/Loss"){
            points(cnvobj$pos[SNPid],rep(major[j],length(SNPid)),type="l",col=gaincol,lwd=linewidth)
            points(rep(cnvobj$pos[mydata$SNP.begin[j]],2),c(mu0,major[j]),type="l",col=gaincol)
            points(rep(cnvobj$pos[mydata$SNP.end[j]],2),c(mu0,major[j]),type="l",col=gaincol)
          }else if (Type[j]=="Loss/Loss"){
            points(cnvobj$pos[SNPid],rep(major[j],length(SNPid)),type="l",col=losscol,lwd=linewidth)
            points(rep(cnvobj$pos[mydata$SNP.begin[j]],2),c(major[j],mu0),type="l",col=losscol)
            points(rep(cnvobj$pos[mydata$SNP.end[j]],2),c(major[j],mu0),type="l",col=losscol)
          }
          if (Type[j]=="Gain/Loss" || Type[j]=="Normal/Loss" || Type[j]=="Loss/Loss"){
            points(cnvobj$pos[SNPid],rep(minor[j],length(SNPid)),type="l",col=losscol,lwd=linewidth)
            points(rep(cnvobj$pos[mydata$SNP.begin[j]],2),c(mu0,minor[j]),type="l",col=losscol)
            points(rep(cnvobj$pos[mydata$SNP.end[j]],2),c(mu0,minor[j]),type="l",col=losscol)
          }else if (Type[j]=="Gain/Gain"){ 
            points(cnvobj$pos[SNPid],rep(minor[j],length(SNPid)),type="l",col=gaincol,lwd=linewidth)
            points(rep(cnvobj$pos[mydata$SNP.begin[j]],2),c(minor[j],mu0),type="l",col=gaincol)
            points(rep(cnvobj$pos[mydata$SNP.end[j]],2),c(minor[j],mu0),type="l",col=gaincol)
          }
        }
      }
    }else{
      plot(x,normal,type="l",col=normalcol,lwd=linewidth,xlab="SNP index",ylab="Copy Number",main=mainlab,ylim=c(0,max(2,max(major))))
      if (n>0){
        for (j in 1:n){
          SNPid = max(a,mydata$SNP.begin[j]):min(b,mydata$SNP.end[j])
          if (Type[j] == "Gain/Gain" || Type[j]=="Gain/Normal" || Type[j]=="Gain/Loss"){
            points(SNPid,rep(major[j],length(SNPid)),type="l",col=gaincol,lwd=linewidth)
            points(rep(mydata$SNP.begin[j],2),c(mu0,major[j]),type="l",col=gaincol)
            points(rep(mydata$SNP.end[j],2),c(mu0,major[j]),type="l",col=gaincol)
          }else if (Type[j]=="Loss/Loss"){
            points(SNPid,rep(major[j],length(SNPid)),type="l",col=losscol,lwd=linewidth)
            points(rep(mydata$SNP.begin[j],2),c(major[j],mu0),type="l",col=losscol)
            points(rep(mydata$SNP.end[j],2),c(major[j],mu0),type="l",col=losscol)
          }
          if (Type[j]=="Gain/Loss" || Type[j]=="Normal/Loss" || Type[j]=="Loss/Loss"){
            points(SNPid,rep(minor[j],length(SNPid)),type="l",col=losscol,lwd=linewidth)
            points(rep(mydata$SNP.begin[j],2),c(mu0,minor[j]),type="l",col=losscol)
            points(rep(mydata$SNP.end[j],2),c(mu0,minor[j]),type="l",col=losscol)
          }else if (Type[j]=="Gain/Gain"){ 
            points(SNPid,rep(minor[j],length(SNPid)),type="l",col=gaincol,lwd=linewidth)
            points(rep(mydata$SNP.begin[j],2),c(minor[j],mu0),type="l",col=gaincol) 
            points(rep(mydata$SNP.end[j],2),c(minor[j],mu0),type="l",col=gaincol)
          }
        }
      }
    } 
}

