`intmode` <-
function(x){
  xt<-tabulate(x)
  xmode<-which(xt == max(xt))
  xmode
}

