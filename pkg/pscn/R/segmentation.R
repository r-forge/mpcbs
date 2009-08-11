`segmentation` <-
function(samplename, chrs=1:23, gender=NULL, MIN.SNPS=20, regroup.percent=-1, combine.alpha=0.05, ...){
    for (i in chrs){
      load(paste(samplename=samplename,".Chr",i,".Smooth.Rdata",sep=""))    
      cnvobj = segment(cnvobj, MIN.SNPS=20, ...)
        cnvobj = suppressWarnings(regroup.theta(cnvobj, level0=level0, level1=level1, percent=regroup.percent))
        if (!is.null(cnvobj$chpts)){
        cnvobj = combine(cnvobj, alpha=combine.alpha)
        if (length(cnvobj$chpts)!=length(cnvobj$chpts.new)){
          n1 = length(cnvobj$chpts.new)
            while(length(cnvobj$chpts.new)>0){
                cnvobj = suppressWarnings(regroup.theta(cnvobj, level0=0.1, level1=0.9, percent=regroup.percent))
                cnvobj = combine(cnvobj, alpha=combine.alpha)
                n2 = length(cnvobj$chpts.new)
                if (n1==n2) break
                n1 = n2
            }
          rm(n1)
        }
      }
    save(cnvobj, file=paste(cnvobj$label,".Segment.Rdata",sep=""))    
  }
  group.number(samplename=samplename, chrs=chrs, gender=gender, MIN.SNPS=MIN.SNPS)
}

