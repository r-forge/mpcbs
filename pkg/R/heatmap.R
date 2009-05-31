`heatmap` <-
function(y,chrom=NULL,loc=NULL, horizontal=TRUE, pos.string="Probe order",...){
    nsamples = ncol(y)
    if(is.null(loc)) loc=1:nrow(y)
    y = y[loc,]
    chrom = chrom[loc]
    
    T = nrow(y)

    
    if(!is.null(chrom)){
      temp=chrom[2:length(chrom)] - chrom[1:(length(chrom)-1)]
      chrom.start = which(temp==1)+1
      final = T
      chromnums = chrom[chrom.start]
      chrom.start = c(0, chrom.start, final)
      chromnums = c(chromnums[1]-1,chromnums)
    } else {
      chrom.start=vector(length=0)
    }
    
    if(horizontal){
      image.plot(loc,1:ncol(y),y,xlab=pos.string,ylab="Sample",...)

      if(length(chrom.start)>1){
        for(i in 1:(length(chrom.start)-1)){
          segments(chrom.start[i],0,chrom.start[i],nsamples)
          mtext(format(chromnums[i]),side=3,at=(chrom.start[i]+chrom.start[i+1])/2)
          if(i==length(chrom.start)-1){
            segments(chrom.start[length(chrom.start)],0,chrom.start[length(chrom.start)],nsamples)
          }
        }
      }
    } else {
      image.plot(1:ncol(y),loc,t(y),ylab=pos.string,xlab="Sample",...)
      if(length(chrom.start)>1){
        for(i in 1:(length(chrom.start)-1)){
          segments(0, chrom.start[i],nsamples, chrom.start[i])
          mtext(format(chromnums[i]),side=4,at=(chrom.start[i]+chrom.start[i+1])/2)
          if(i==length(chrom.start)-1){
            segments(0,chrom.start[length(chrom.start)],nsamples,chrom.start[length(chrom.start)])
          }
        }
      }
    }
}

