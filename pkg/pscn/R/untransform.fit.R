`untransform.fit` <-
function(cnvobj, seg=TRUE){
                                        # Assumes that the X,B values were obtained from logR, theta values from getXY, and
                                        # NO FURTHER PROCESSING WAS DONE.

  if(seg && !is.null(cnvobj$sig.seg)){
    X = cnvobj$sig.seg[,1]
    Y = cnvobj$sig.seg[,2]
    alpha = cnvobj$theta.seg[,1]
    beta = cnvobj$theta.seg[,2]
  } else {
    X =  cnvobj$sig[,1]
    Y = cnvobj$sig[,2]
    alpha = cnvobj$theta[,1]
    beta = cnvobj$theta[,2]    
  }
  temp1=illumina.getlogRbfreq.fromXY(X,Y,cnvobj$illumina$OFFSET)
  logR.fit = temp1$logR
  bfreq.fit = temp1$bfreq

  cnvobj$illumina$logRfit = logR.fit
  cnvobj$illumina$bfreqfit = bfreq.fit

  major.tr = pmax(alpha+beta, alpha-beta)
  minor.tr = pmin(alpha+beta, alpha-beta)
  temp1=illumina.getlogRbfreq.fromXY(major.tr,minor.tr,cnvobj$illumina$OFFSET)
  temp2=illumina.getAB(temp1$logR, temp1$bfreq)
  cnvobj$illumina$majorfit = temp2$B
  cnvobj$illumina$minorfit = temp2$A

  temp3 = illumina.getAB(logR.fit, bfreq.fit)
  cnvobj$illumina$Afit = temp3$A
  cnvobj$illumina$Bfit = temp3$B
  cnvobj
}

