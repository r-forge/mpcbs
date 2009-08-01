`centro` <-
function(cnvobj, MIN.GAP = 2e6){
  pos = cnvobj$pos
  N = length(pos)
  Pos.begin = pos[1]
  Pos.end = pos[N]
  Effective.pos = Pos.end - Pos.begin
  Centromere = 0
  pos.diff = pos[2:N]-pos[1:(N-1)]
  max.pos.diff = max(pos.diff)
  SNP.star = NA
  if (max.pos.diff > MIN.GAP){
    SNP.star = which(pos.diff==max.pos.diff)
    Centromere = paste(pos[SNP.star],pos[SNP.star+1],sep="-")
    Effective.pos = Effective.pos - max.pos.diff
  }
  return(c(Pos.begin, Pos.end, Effective.pos, SNP.star))
}

