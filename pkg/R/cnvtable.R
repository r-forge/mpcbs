`cnvtable` <-
function(chrom,pos,featuremat,segments,sampleswithsegment){
    nsegs = nrow(segments)
    
    NCOLS = 11
    cnvlist=matrix(0,nrow=0,ncol=NCOLS)
   
    for(i in 1:nsegs){
        st = segments[i,1]; ed = segments[i,2]
        carriers = which(sampleswithsegment[,i]==1)
        n.with.cnv = length(carriers)
        carriers.gain = carriers[which(featuremat[carriers,i]>0)]
        carriers.loss = carriers[which(featuremat[carriers,i]<0)]
        if(length(carriers.gain)>0) mean.gain = mean(featuremat[carriers.gain,i])
        else mean.gain=0
        if(length(carriers.loss)>0) mean.loss = mean(featuremat[carriers.loss,i])
        else mean.loss=0
        cnvlist = rbind(cnvlist,c(st,ed,ed-st+1,chrom[st],pos[st],pos[ed],length(carriers),length(carriers.gain),mean.gain,length(carriers.loss),mean.loss))
    }
    cnvlist=as.data.frame(cnvlist)
    names(cnvlist)<-c("Start index","End index","# SNPs","Chromosome","Start pos","End pos","# Carriers","# Gain","Mean gain","# Loss","Mean loss")
  
    cnvlist
}

