`print.cnv2d` <-
function(cnvobj){
    cat(paste("\n\nCNV object ",cnvobj$label,": total ", length(cnvobj), " SNPs on ",cnvobj$platform," platform.\n",sep=""))
    if(!is.null(cnvobj$sig)) cat("BCMIX fitted values available.\n")
    if(!is.null(cnvobj$sig.seg)) cat(paste("Segmented: total ",length(cnvobj$chpts)," change-points.\n\n",sep="")) 
    
    if(!is.null(cnvobj$cnvtable))  print(format(cnvobj$cnvtable,digits=2))        
}

