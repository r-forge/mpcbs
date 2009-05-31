`getCutoffMultisampleWeightedChisq` <-
function(pval,m,delta,win,N,alpha){
    
    cat("Computing threshold for weighted chi-square...\n")
    THRES = 0.1*pval
    currb = 1
    prevsmallerb = currb
    prevlargerb = 200
    currpval=1
    
    while( abs(currpval-pval)>THRES ){
#        cat("pval =", pval, ", currpval = ",currpval,", THRES=", THRES,".\n",sep="")
        
        if( currpval>pval){
            # need to increase b.
            prevsmallerb = currb
            currb = currb+(prevlargerb-currb)/2
        } else {
            # need to decrease b.
            prevlargerb = currb
            currb = currb - (currb-prevsmallerb)/2
        }
    
#        cat("currb = ",currb,"\n")
        currpval = pvalueMultisampleWeightedChisq(currb,m,delta,win,alpha,N);
    }    
    currb
}

