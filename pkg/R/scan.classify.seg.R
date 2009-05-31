`scan.classify.seg` <-
function(this.y,seg,y.var,CHISQ.PVAL.THRESH=0.001, MIN.REQ.ABSDIFF=0, MIN.SUFF.ABSDIFF=Inf){
        
        sample.pval = computeZ.squarewave.sample(t(this.y),seg,y.var)$pval
        if(seg[2]-seg[1]>1){ 
            sample.med.diff=apply(as.matrix(this.y[(seg[1]+1):seg[2],]),2,median) - apply(as.matrix(this.y[-c((seg[1]+1):seg[2]),]),2,median)
        } else { 
            sample.med.diff=as.matrix(this.y[(seg[1]+1):seg[2],]) - apply(as.matrix(this.y[-c((seg[1]+1):seg[2]),]),2,median)
        }
       
        sample.abs.med.diff = abs(sample.med.diff)
        pass1 = sample.pval<CHISQ.PVAL.THRESH
        pass2 = sample.abs.med.diff > MIN.SUFF.ABSDIFF
        cut1 = sample.abs.med.diff  < MIN.REQ.ABSDIFF

        sampleswithsegment = (pass1 | pass2) & (!cut1)
        sampleswithsegment

}

