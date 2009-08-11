`plot.sig` <-
function(x, loc=NULL, which.chrom=NA, which.plot, seg, plot.pos=FALSE, display.legend=TRUE, main=NULL,color.points=FALSE,...){
   if (!inherits(x, 'cnv2d')) 
      stop("First argument must be a cnv2d object.")
  
    n=dim(x$intensity)[1]

   if(is.null(main)){
     main = x$samplename
   }


   if(color.points){
     raw.colors=c("pink","skyblue","green","gray") 
   } else {
     raw.colors=rep("gray",4)
   }
   
    if (which.plot=="sig"){
        plotallele1= TRUE
        plotallele2=TRUE
    }
    
    if(which.plot =="allele1"){
        plotallele1 = TRUE
        plotallele2 = FALSE
    }
    if(which.plot =="allele2"){
        plotallele1 = FALSE
        plotallele2 = TRUE
    }

    if(seg && !is.null(x$sig.seg)){
        sig = x$sig.seg
    } else {
        sig =x$sig
    }

    if(plot.pos){
        xval = x$pos
        xlabel = "Position (bp)"
    } else {
        xval = c(1:length(x))
        xlabel = "SNP index"
    }

   if(!is.na(which.chrom)) {
     if(is.null(x$mixstate)) mixstateloc = rep(1,length(loc))
     else mixstateloc = x$mixstate[ind]
     ind = which(x$chromosome == which.chrom)

        if (plotallele1){
            plot(xval[ind], x$intensity[ind,1], xlab=xlabel,ylab="B'", pch=4,cex=1.5,
                col=raw.colors[mixstateloc], main=main,...)

            if(!is.null(sig)){
                points(xval[which(x$mixstate==1  & x$chromosome == which.chrom)],sig[which(x$mixstate==1  & x$chromosome == which.chrom),1],col="red")
                points(xval[which(x$mixstate==2  & x$chromosome == which.chrom)],sig[which(x$mixstate==2  & x$chromosome == which.chrom),1],col="blue")
                points(xval[which(x$mixstate==3  & x$chromosome == which.chrom)],sig[which(x$mixstate==3  & x$chromosome == which.chrom),1],col="darkgreen")
                points(xval[which(x$mixstate==4  & x$chromosome == which.chrom)],sig[which(x$mixstate==4  & x$chromosome == which.chrom),1],col="black")
            }
        }

        if(plotallele2){
            plot(xval[ind], x$intensity[ind,2], xlab=xlabel,ylab="A'", pch=4, cex=1.5,
                col=raw.colors[mixstateloc],...)


            if(!is.null(sig)){
                points(xval[which(x$mixstate==1  & x$chromosome == which.chrom)],sig[which(x$mixstate==1  & x$chromosome == which.chrom),2],col="red")
                points(xval[which(x$mixstate==2  & x$chromosome == which.chrom)],sig[which(x$mixstate==2  & x$chromosome == which.chrom),2],col="blue")
                points(xval[which(x$mixstate==3  & x$chromosome == which.chrom)],sig[which(x$mixstate==3  & x$chromosome == which.chrom),2],col="darkgreen")
                points(xval[which(x$mixstate==4  & x$chromosome == which.chrom)],sig[which(x$mixstate==4  & x$chromosome == which.chrom),2],col="black")
            }
        }
   } else { 
   
        if(is.null(loc)){
            loc = c(1:length(x))
        }
        
        
        
        if(plotallele1){

          sigloc = sig[loc,]
          if(is.null(x$mixstate)) mixstateloc = rep(1,length(loc))
          else mixstateloc = x$mixstate[loc]
          plot(xval[loc], x$intensity[loc,1], xlab=xlabel,ylab="B'",  pch=4, cex=1.5,
                col=raw.colors[mixstateloc], main=main,...)

            
            if(!is.null(sig)){
                points(xval[loc[which(mixstateloc==1)]],sigloc[which(mixstateloc==1),1],col="red")
                points(xval[loc[which(mixstateloc==2)]],sigloc[which(mixstateloc==2),1],col="blue")
                points(xval[loc[which(mixstateloc==3)]],sigloc[which(mixstateloc==3),1],col="darkgreen")
                points(xval[loc[which(mixstateloc==4)]],sigloc[which(mixstateloc==4),1],col="black")
            }

            maxy = max(x$intensity[loc,1])
            miny = min(x$intensity[loc,1])
            for(i in loc[2:length(loc)]){
                if(x$chromosome[i]!=x$chromosome[i-1]){
                    segments(xval[i]-0.5, maxy, xval[i]-0.5, miny)
                }
            }

        }
        
        
        if(plotallele2){
            sigloc = sig[loc,]
            if(is.null(x$mixstate)) mixstateloc = rep(1,length(loc))
            else mixstateloc = x$mixstate[loc]

          plot(xval[loc], x$intensity[loc,2], xlab=xlabel,ylab="A'",  pch=4, cex=1.5,
                col=raw.colors[mixstateloc], main=main,...)

            
            if(!is.null(sig)){
                points(xval[loc[which(mixstateloc==1)]],sigloc[which(mixstateloc==1),2],col="red")
                points(xval[loc[which(mixstateloc==2)]],sigloc[which(mixstateloc==2),2],col="blue")
                points(xval[loc[which(mixstateloc==3)]],sigloc[which(mixstateloc==3),2],col="darkgreen")
                points(xval[loc[which(mixstateloc==4)]],sigloc[which(mixstateloc==4),2],col="black")
            }

            maxy = max(x$intensity[loc,2])
            miny = min(x$intensity[loc,2])
            for(i in loc[2:length(loc)]){
                if(x$chromosome[i]!=x$chromosome[i-1]){
                    segments(xval[i]-0.5, maxy, xval[i]-0.5, miny)
                }
            }
        }
   }
}

