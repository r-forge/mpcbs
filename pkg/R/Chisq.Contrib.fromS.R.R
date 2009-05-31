`Chisq.Contrib.fromS.R` <-
function(S,SST,imap,chpts){
    T = nrow(S)
    N = ncol(S)
    dfnum <- 1;
    n.segs <- nrow(chpts)
    
    ret.m <- matrix(nrow=n.segs, ncol=N, data=0)
 
    totalsnps = imap[T,]-imap[1,]
    
    for(i in 1:n.segs){
        st <- chpts[i,1]
        ed <- chpts[i,2]
        w=ed-st
        
        for(j in 1:N){
            nsnps = imap[ed,j]-imap[st,j]        
            mn1 <- S[T,j]/totalsnps[j]
            SSb <- (S[ed,j]-S[st,j] -nsnps*mn1)^2 / (nsnps*(1-nsnps/totalsnps[j]));
            SSw = SST[j]-SSb;
            ret.m[i,j] = (SSb/dfnum)/(SSw/(totalsnps[j]-2));
        }

    }
    
    return(ret.m)
}

