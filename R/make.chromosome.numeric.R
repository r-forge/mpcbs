`make.chromosome.numeric` <-
function(chrom){

  if(!is.numeric(chrom)){
    chromosome=rep(0,length(chrom))
    chromosome[chrom=="1"]=1
    chromosome[chrom=="2"]=2
    chromosome[chrom=="3"]=3
    chromosome[chrom=="4"]=4
    chromosome[chrom=="5"]=5
    chromosome[chrom=="6"]=6
    chromosome[chrom=="7"]=7
    chromosome[chrom=="8"]=8
    chromosome[chrom=="9"]=9
    chromosome[chrom=="10"]=10
    chromosome[chrom=="11"]=11
    chromosome[chrom=="12"]=12
    chromosome[chrom=="13"]=13
    chromosome[chrom=="14"]=14
    chromosome[chrom=="15"]=15
    chromosome[chrom=="16"]=16
    chromosome[chrom=="17"]=17
    chromosome[chrom=="18"]=18
    chromosome[chrom=="19"]=19
    chromosome[chrom=="20"]=20
    chromosome[chrom=="21"]=21
    chromosome[chrom=="22"]=22
    chromosome[chrom=="X"]=23
    chromosome[chrom=="XY"]=24
    chromosome[chrom=="Y"]=25  
  }
  chromosome
}

