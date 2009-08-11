`NA.summary` <-
function(mat){
  natype = c("NA.0.1","NA.0.0","NA.-1.1","NA.-1.0")
  leng = rep(0,4)
  idtotal = c()
  for (i in 1:4){
    id = which(mat$Type==natype[i])
    leng[i] = sum(length(id))
    idtotal = c(idtotal,id)
  }
  list(leng=leng,idtotal=idtotal)
}

