`illumina.getAB` <-
function(logR, bfreq, R.BASELINE=2){
  # IMPORTANT: assumes that logR=0 corresponds to 2 copies.
  # That is, R=exp(R)=1 corresponds to 2 copies.  If not,
  # need to alter R.BASELINE.
  R = exp(logR)
    R = R*R.BASELINE
    A=R/(1+tan(bfreq*pi/2))
    B=R-A
    list(A=A,B=B)
}

