`create.cnv2d` <-
function(obs, samplename, chromosome=NULL, pos=NULL, label=NULL, platform="Unknown")
{        
    dimobs=dim(obs)
    n=dimobs[1]

    if(is.null(chromosome)) chromosome = rep(1,n)                                     
    if(is.null(pos)) pos = c(1:n)

    if (n != length(chromosome) || n != length(pos))
        stop("intensity, chromosome, and pos need to have same length.")

    isvalid = which(!(is.na(obs[,1]) | is.na(obs[,2])))    
    
    n=length(isvalid)
    nchroms = length(unique(chromosome[isvalid]))

    if(platform=="Illumina" || platform=="Affymetrix" || platform=="Unknown"){
        logR = obs[isvalid,1]
        theta = obs[isvalid,2]

        if(max(sum(is.na(logR)), sum(is.na(theta))) > 0){
          cat("There are ",sum(is.na(logR)), " missing log R values and ", sum(is.na(theta)), " missing b freq values.  Replace by NORMAL values.")
          logR[is.na(logR)] =0
          theta[is.na(logR)] =0
        }

        temp = illumina.getXY(logR,theta)
        OFFSET = temp$OFFSET
        AB = illumina.getAB(logR,theta)
        A=AB$A
        B=AB$B
        illumina = list(logR = logR,theta=theta,OFFSET=OFFSET, A=A, B=B) # The R, theta values must be what was originally passed in.
        intensity=cbind(temp$X, temp$Y)
      } else {                                                                                                  
        intensity=obs[isvalid,]
        illumina = NULL
    }

    if(is.null(label)){
      if(is.null(chromosome)){
        label= samplename
      }else{
        label = paste(samplename,".Chr",chromosome[1],sep="")
      }
    }
    
    cnvobj <-
        list(intensity = intensity,
             chromosome = chromosome[isvalid],
             pos= pos[isvalid],
             nchroms = nchroms,
             samplename = samplename,
             label = label,
             statemat = vector("list",4),
             illumina = illumina,
             platform = platform  
            )
    cnvobj$statemat[[1]] = matrix(c(2,0,0,0),2,2,byrow=TRUE)
    cnvobj$statemat[[2]] = matrix(c(1,1,1,-1),2,2,byrow=TRUE)
    cnvobj$statemat[[3]] = matrix(c(1,1,-1,1),2,2,byrow=TRUE)
    cnvobj$statemat[[4]] = matrix(c(0,2,0,0),2,2,byrow=TRUE)
 
    class(cnvobj) <- "cnv2d"
    attr(cnvobj, "call") <- sys.call()    
    inherits(cnvobj,"cnv")
    
    cnvobj$intensity.unnormed=cnvobj$intensity
    cnvobj$normparams$scale.fac=c(1,1)
    cnvobj$normparams$cstar=c(1,1)

    cat(paste("CNV object ",label," created: Number of probes:", n, ", number of chromosomes:", nchroms, ".\n",sep=""))
    
    cnvobj
}

