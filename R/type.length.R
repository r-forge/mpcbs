`type.length` <-
function(mat){
  Type = c("Gain/Gain","1.Gain/Gain","t.Gain/Gain", "Gain/Normal","t.Gain/Normal","Gain/Loss","t.Gain/Loss","Normal/Loss","t.Normal/Loss","Loss/Loss","1.Loss/Loss","t.Loss/Loss")
  leng = rep(0,12)
  for (i in 1:12){
    id = which(mat$Type==Type[i])
    temp = as.numeric(mat$Pos.length[id])
    leng[i] = sum(temp)
  }
  leng2 = rep(0,5)
  leng2[1] = sum(leng[1:3])
  leng2[2] = sum(leng[4:5])
  leng2[3] = sum(leng[6:7])
  leng2[4] = sum(leng[8:9])
  leng2[5] = sum(leng[10:12])
  leng2
}

