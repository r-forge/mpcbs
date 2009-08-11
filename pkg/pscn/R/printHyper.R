`printHyper` <-
function(hyper,filename=NULL){
    if(!is.null(filename)) sink(filename)
    cat(format(hyper$p,scientific=FALSE), " ",sep="")
    cat(format(hyper$a,scientific=FALSE), " ",sep="")
    cat(format(hyper$b,scientific=FALSE), "\n",sep="")
    cat(hyper$mu[1]," ",hyper$mu[2],"\n",sep="")
    cat(hyper$v[1,1]," ",hyper$v[1,2],"\n",hyper$v[2,1]," ",hyper$v[2,2],"\n",sep="")
    cat(hyper$sigAA[1,1]," ",hyper$sigAA[1,2],"\n",hyper$sigAA[2,1]," ",hyper$sigAA[2,2],"\n",sep="")
    cat(hyper$sigAB[1,1]," ",hyper$sigAB[1,2],"\n",hyper$sigAB[2,1]," ",hyper$sigAB[2,2],"\n",sep="")
    cat(hyper$sigBA[1,1]," ",hyper$sigBA[1,2],"\n",hyper$sigBA[2,1]," ",hyper$sigBA[2,2],"\n",sep="")
    cat(hyper$sigBB[1,1]," ",hyper$sigBB[1,2],"\n",hyper$sigBB[2,1]," ",hyper$sigBB[2,2],"\n",sep="")
    cat(hyper$base[1]," ",hyper$base[2],"\n",sep="")

    if(!is.null(filename)) sink()
}

