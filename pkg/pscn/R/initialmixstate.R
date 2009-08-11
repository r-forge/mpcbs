`initialmixstate` <-
function(cnvobj){
  n = length(cnvobj$illumina$theta)
  mixstate = rep(1,n)
  id2 = which(cnvobj$illumina$theta<0.99)
  mixstate[id2] = 2
  id3 = which(cnvobj$illumina$theta<0.5)
  mixstate[id3] = 3
  id4 = which(cnvobj$illumina$theta<0.01)
  mixstate[id4] = 4
  mixstate
}

