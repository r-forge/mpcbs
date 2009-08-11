`segment` <-
function(cnvobj, verbose=FALSE, ...){

    if(verbose) cat("\nSegmenting ",cnvobj$label,"...\n",sep="")
    
    cnvobj=segment.sub(cnvobj, verbose, start.from.seg=FALSE,...)
 
    while(TRUE){
        if(verbose) cat("Round 1 error checking...")
        old.cp = cnvobj$chpts
        cnvobj = segment.sub(cnvobj, verbose, start.from.seg=TRUE,...)
        new.cp = cnvobj$chpts
        if(length(old.cp)-length(new.cp)==0) break        
 # Hao 09/09/2008
        if(is.null(cnvobj$chpts)) break
    }
       
   if(verbose) cat("Done.\n\n")
    cnvobj
}

