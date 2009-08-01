`illumina.getlogRbfreq` <-
function(A,B,OFFSET=2){
    logR = log((A+B)/OFFSET)
    bfreq = atan(B/A)/(pi/2)
    list(logR=logR,bfreq=bfreq)
}

