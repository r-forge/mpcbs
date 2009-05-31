`Chisq.Contrib.R` <-
function(y, T, chpts){
    dfnum <- 1;
    dfden <- T-2;

    N <- nrow(y)
    ret.m <- matrix(nrow=nrow(chpts), ncol=N, data=0)
    # chisq holds the contribution of each sample to the aberration.

    n.segs <- nrow(chpts)

    for(i in 1:n.segs){
        st <- chpts[i,1];
        ed <- chpts[i,2];
        w <- ed-st; # Changed 11/4, previously w<-ed-st+1.

        for(j in 1:N){
            mn1 <- sum(y[j,])/T
            SST <- sum( (y[j,]- mn1)^2 );
            SSb <- (( sum(y[j,(st+1):ed])-w*mn1)^2) / (w*(1-w/T));
            SSw = SST-SSb;
            ret.m[i,j] = (SSb/dfnum)/(SSw/dfden);
        }

    }
    return(ret.m)
}

