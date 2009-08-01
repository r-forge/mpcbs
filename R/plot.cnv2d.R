`plot.cnv2d` <-
function (x, loc=NULL, which.chrom=NA, which.plot = "sig", 
                        seg = TRUE, display.legend=TRUE,color.points=FALSE, ...) {
    
    if(which.plot=="sig" || which.plot=="allele1" || which.plot=="allele2" ){
        plot.sig(x,loc=loc, which.chrom=which.chrom,which.plot=which.plot, seg=seg, display.legend=display.legend, color.points=color.points,...)
    } else {
       
        plot.onefeature(x,loc=loc, which.chrom=which.chrom, which.plot=which.plot, seg=seg,  display.legend=display.legend,...)
    }
}

