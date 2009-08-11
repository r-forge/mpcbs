`smooth.cnv2d` <-
function(cnvobj,hyper=NULL,mixstate=NULL, niters=1,nBCMIXiters=10,K=20,M=10,outerloop=2, selectHyper=4, verbose=TRUE, plots=FALSE, save.interim=FALSE){

  # UNCOMMENT THESE IF NOT DEBUGGING!!!
#  hyper=NULL; mixstate=NULL; niters=1; nBCMIXiters=5; K=20; M=10; outerloop=2; selectHyper=4; verbose=TRUE; plots=FALSE

  obs = cnvobj$intensity

    # ---- STEP 1: Infer initial state assignment.

    if(is.null(mixstate)){
        mixstate = initialmixstate(cnvobj)
    }
    
    # ---- STEP 2: Set hyper.
    if(!is.null(hyper)) {
        cnvobj$hyper <- hyper
    } else {
        cnvobj$hyper <- initialHyper.cnv2d(cnvobj)
    }

    # ---- STEP 3: Call cppEstHyper.  
    #       Each iteration of cppEstHyper does "niters" iterations of BCMIX/FindState/Hyperparameter estimation.
    #       The number of BCMIX/FindState iterations within each iteration of the above is nBCMIXiters.
    #       First, everything needs to be vectorized for passing to C++.
    
    obsvec = rep(0,nrow(obs)*2)
    obsvec[seq(1,nrow(obs)*2-1,2)] = obs[,1]
    obsvec[seq(2,nrow(obs)*2,2)] = obs[,2]

    cnvobj$fitparams = list(niters=niters,nBCMIXiters=nBCMIXiters,K=K,M=M,outerloop=outerloop, selectHyper=selectHyper)

    for(iter in c(1:outerloop)){
       
       if(verbose) cat("Iteration ",iter,"...")
        # For debugging purpose:
#        if(verbose){
#            cat("\n\nIteration ",iter,", current hyperparameters are:\n")
#            printHyper(cnvobj$hyper)
#        }

        # Call cppEstHyper
        fit =.C("cppEstHyper",   as.double(obsvec), mixstate=as.integer(mixstate), as.integer(nrow(obs)),
                                    newp=as.double(cnvobj$hyper$p), newa=as.double(cnvobj$hyper$a), newb=as.double(cnvobj$hyper$b), 
                                    newbase1=as.double(cnvobj$hyper$base[1]), newbase2=as.double(cnvobj$hyper$base[2]),
                                    newmuvec=as.double(cnvobj$hyper$mu), newvvec=as.double(cnvobj$hyper$vvec),
                                    newsigAAvec=as.double(cnvobj$hyper$sigAAvec), newsigABvec=as.double(cnvobj$hyper$sigABvec),
                                    newsigBAvec=as.double(cnvobj$hyper$sigBAvec), newsigBBvec=as.double(cnvobj$hyper$sigBBvec), 
                                    as.integer(K), as.integer(M), as.integer(niters), as.integer(nBCMIXiters),
                                    estSig=double(length(obsvec)), estBS=double(nrow(obs)), as.integer(selectHyper),
                                    PACKAGE="pscn")     
    
        
        # ---- Unpack everything.
        if(verbose) cat("done.\n")
    
        theta = matrix(fit$estSig,nrow=nrow(obs),ncol=2,byrow=TRUE)
        mixstate=fit$mixstate
        cnvobj$hyper = createHyper(fit$newp, fit$newa, fit$newb, fit$newmuvec, fit$newvvec, fit$newsigAAvec, fit$newsigABvec, fit$newsigBAvec, fit$newsigBBvec, c(fit$newbase1,fit$newbase2))
        sig = matrix(0,nrow(obs),2)
        sig[mixstate==1,] = theta[mixstate==1,]%*%cnvobj$statemat[[1]]
        sig[mixstate==2,] = theta[mixstate==2,]%*%cnvobj$statemat[[2]]
        sig[mixstate==3,] = theta[mixstate==3,]%*%cnvobj$statemat[[3]]
        sig[mixstate==4,] = theta[mixstate==4,]%*%cnvobj$statemat[[4]]
        cnvobj$sig = sig
        cnvobj$theta = theta
        cnvobj$mixstate = mixstate
        cnvobj$estBS = fit$estBS
        
        
        # ---- Save and do diagnostic plots.
        if (iter<outerloop){
          if (save.interim){      
            if(verbose) cat("Saving current state....\n")
            save(cnvobj, file=paste(cnvobj$label,".fit",iter,".Rdata",sep=""))
          }
        }else{
          if(verbose) cat("Saving current state....\n")
          save(cnvobj, file=paste(cnvobj$label,".Smooth.Rdata",sep=""))
        }
        
        if(plots){
            cat("Drawing diagnostic plot....\n")
#            png(paste(cnvobj$label,iter,"progress.png",sep="_"), height=2000, width=2000)
            png(paste(cnvobj$label,iter,"progress.png",sep="_"))
            par(mfrow=c(2,1))
            plot(cnvobj,color.points=TRUE)
            dev.off()
        }
        
    }
    cnvobj
}

