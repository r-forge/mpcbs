`smoothing` <-
function(samplename, inputdata, platform="Unknown", sm.plot=FALSE, ...){
  tempchr = inputdata$Chr
  for (j in 1:23){
    id = length(which(tempchr==j))
    if (id>0){
      chrdata = inputdata[which(tempchr==j),]
      chrdata = chrdata[sort(chrdata$Position,index.return=TRUE)$ix,]
      cnvobj = create.cnv2d(chrdata[,3:4], chromosome=chrdata$Chr, pos=chrdata$Position, samplename=samplename, platform=platform)
      cnvobj = smooth.cnv2d(cnvobj, plots=sm.plot, ...)
    }
    rm(id)
  }
}

