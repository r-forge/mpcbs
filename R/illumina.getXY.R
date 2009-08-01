`illumina.getXY` <-
function(logR,bfreq, OFFSET=-min(logR)){
# The larger the offset, the more weight is given to variation in bfreq and the less is given to variation in logR.
# Since logR is more noisy than bfreq, smaller offsets make the data look more noisy.
# But then again, smaller offsets look more normal.

    R2=(logR+OFFSET)/OFFSET
    X=R2*bfreq
    Y=R2-X
    list(X=X,Y=Y,OFFSET=OFFSET)
}

