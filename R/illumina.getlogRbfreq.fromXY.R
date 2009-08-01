`illumina.getlogRbfreq.fromXY` <-
function(X,Y,OFFSET){
  logR = X+Y-OFFSET
  bfreq = X/(logR+OFFSET)
  list(logR=logR, bfreq=bfreq)
}

