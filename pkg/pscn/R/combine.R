`combine` <-
function(cnvobj, alpha=0.05){
  cnvobj = combine.sub(cnvobj)
  if (!is.null(cnvobj$chpts.new)){
    n1 = length(cnvobj$chpts.new)
    while(length(cnvobj$chpts.new)>0){
      cnvobj = combine.sub(cnvobj, alpha=alpha)
      n2 = length(cnvobj$chpts.new)
      if(n1==n2){break}
      n1 = n2
    }
  }
  cnvobj
}

