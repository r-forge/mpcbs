`cnvlist.update` <-
function(mydata, alpha=0.01, threshold = 0.18, hard.threshold=0.1){
  n = dim(mydata)[1]
  p = alpha/(2*n)
  mu0 = mydata$Normal.copy[1]
  mydata$Type = as.character(mydata$Type)
  for (j in 1:n){
    Mp = mydata$major.pvalue[j]
    mp = mydata$minor.pvalue[j]
    M = as.numeric(strsplit(as.character(mydata$Value[j]),"/")[[1]][1])
    m = as.numeric(strsplit(as.character(mydata$Value[j]),"/")[[1]][2])
    if (mp>(1-p) && (m-mu0)>threshold){
      mydata$Type[j] = "Gain/Gain"
    }else if (Mp<p && (mu0-M)>threshold){
      mydata$Type[j] = "Loss/Loss"
    }else if (Mp>(1-p) && mp<p && (M-mu0)>threshold && (mu0-m)>threshold){
      mydata$Type[j] = "Gain/Loss"
    }else if (Mp>(1-p) && ( (mp>(1-p) && (m-mu0)<threshold) || (mp>p && mp<(1-p)) || (mp<p && (mu0-m)<threshold) && (mu0-m)<(M-mu0) )){
      mydata$Type[j] = "Gain/Normal"
    }else if (mp<p && ( (Mp>(1-p) && (M-mu0)<(mu0-m) && (M-mu0)<threshold) || (Mp>p && Mp<(1-p)) || (Mp<p && (mu0-M)<threshold) )){
      mydata$Type[j] = "Normal/Loss"
    }else{
      mydata$Type[j] = "NA"
    }
  }
  # use hard threshold
  for (j in 1:n){
      M = as.numeric(strsplit(as.character(mydata$Value[j]), 
            "/")[[1]][1])
      m = as.numeric(strsplit(as.character(mydata$Value[j]), 
            "/")[[1]][2])
      if ((M-mu0)<hard.threshold && (mu0-m)<hard.threshold){
        mydata$Type[j] = "NA"
      }
  }

  return(mydata)
}

