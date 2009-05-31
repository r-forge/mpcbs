\name{msscan}
\alias{msscan}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ ~~function to do ... ~~ }
\description{
  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
msscan(y, win, ALPHA = 0, GLOBAL.PVAL.CUTOFF = 1e-04, SAMPLE.PVAL.CUTOFF = 1e-04, SAMPLE.LOW.THRESH = 0.1, SAMPLE.HIGH.THRESH = 0.5, SAMPLE.CUTOFF.METHOD = 4, verbose = TRUE, WCHISQ.CUTOFF = NA, SAMPLE.THRESHOLD = NA, MIN.SNPS = 2, Z = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{y}{ ~~Describe \code{y} here~~ }
  \item{win}{ ~~Describe \code{win} here~~ }
  \item{ALPHA}{ ~~Describe \code{ALPHA} here~~ }
  \item{GLOBAL.PVAL.CUTOFF}{ ~~Describe \code{GLOBAL.PVAL.CUTOFF} here~~ }
  \item{SAMPLE.PVAL.CUTOFF}{ ~~Describe \code{SAMPLE.PVAL.CUTOFF} here~~ }
  \item{SAMPLE.LOW.THRESH}{ ~~Describe \code{SAMPLE.LOW.THRESH} here~~ }
  \item{SAMPLE.HIGH.THRESH}{ ~~Describe \code{SAMPLE.HIGH.THRESH} here~~ }
  \item{SAMPLE.CUTOFF.METHOD}{ ~~Describe \code{SAMPLE.CUTOFF.METHOD} here~~ }
  \item{verbose}{ ~~Describe \code{verbose} here~~ }
  \item{WCHISQ.CUTOFF}{ ~~Describe \code{WCHISQ.CUTOFF} here~~ }
  \item{SAMPLE.THRESHOLD}{ ~~Describe \code{SAMPLE.THRESHOLD} here~~ }
  \item{MIN.SNPS}{ ~~Describe \code{MIN.SNPS} here~~ }
  \item{Z}{ ~~Describe \code{Z} here~~ }
}
\details{
  ~~ If necessary, more details than the description above ~~
}
\value{
  ~Describe the value returned
  If it is a LIST, use
  \item{comp1 }{Description of 'comp1'}
  \item{comp2 }{Description of 'comp2'}
  ...
}
\references{ ~put references to the literature/web site here ~ }
\author{ ~~who you are~~ }
\note{ ~~further notes~~ 

 ~Make other sections like Warning with \section{Warning }{....} ~
}
\seealso{ ~~objects to See Also as \code{\link{help}}, ~~~ }
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function(y,win,ALPHA=0,GLOBAL.PVAL.CUTOFF=0.0001,SAMPLE.PVAL.CUTOFF=0.0001,
                    SAMPLE.LOW.THRESH = 0.1, SAMPLE.HIGH.THRESH=0.5,
                    SAMPLE.CUTOFF.METHOD=4, verbose=TRUE,
                    WCHISQ.CUTOFF=NA, SAMPLE.THRESHOLD=NA,
                    MIN.SNPS = 2, Z=NULL){
    
    N=dim(y)[1]
    T=dim(y)[2]
    DELTA = 1/T
    
    if(is.na(WCHISQ.CUTOFF)){
        WCHISQ.CUTOFF = getCutoffMultisampleWeightedChisq(GLOBAL.PVAL.CUTOFF,T,DELTA,win,N,ALPHA)
    }
    
    if(SAMPLE.CUTOFF.METHOD==2 || SAMPLE.CUTOFF.METHOD==3 || SAMPLE.CUTOFF.METHOD==4){
        if(is.na(SAMPLE.THRESHOLD)){
            SAMPLE.THRESHOLD = getCutoffEpidemicChangeLocationKnown(SAMPLE.PVAL.CUTOFF)
        }
    }
    
    if(verbose){
        cat("Running msscan: ",N," samples, ",T," SNPs, ALPHA=",ALPHA,"\n", sep="")
        cat("Global p-value cutoff: ",GLOBAL.PVAL.CUTOFF,", corresponding Chisq cutoff: ",WCHISQ.CUTOFF, ".\n", sep="");
    }
    
#    if(msscan.debug.trace){
#         print("about to call ComputeZ.c")
#         print(paste("dim(y) = ", nrow(y), ncol(y)) )
#         print(paste("ALPHA = ", ALPHA) )
#         scan()
#    }

    if(is.null(Z)){
        Z = ComputeZ.C(y, T, win, ALPHA)
        Z[,1:MIN.SNPS-1] = 0
#        if(msscan.debug.trace){
#            print("about to call computeMoments")
#            scan()
#        }
        # use this gnorm to compute moments, which assumes u is a gaussian.
        g2<-function(u){ u^2*exp(u^2/2)/(ALPHA+exp(u^2/2))}
        g2.moments=computeMoments(g2,0);
        Z = (Z - N*g2.moments$psidot)/sqrt(g2.moments$psidotdot*N)    
    }
       
    Z.passinds = which(Z>WCHISQ.CUTOFF,arr.ind=TRUE)

    if(length(Z.passinds) == 0){
         y.hat <- y
         y.hat[,] <- mean(y)
         return( list(yhat=yhat, chpts=matrix(nrow=0, ncol=2), sampleswithsegment=NULL,chisq=NULL,Z=NULL)  )
    }

    Z.pass = Z[Z.passinds]
    Z.pass.order =order(Z.pass,decreasing=TRUE)
    Z.pass = Z.pass[Z.pass.order]
    Z.passinds = Z.passinds[Z.pass.order,]
    
    chpts = cbind(Z.passinds[,1]+1,Z.passinds[,1]+Z.passinds[,2])
    chpts.keep = Remove.Overlap.C(chpts)
    chpts = chpts[chpts.keep,]    

    if(length(chpts)==2) chpts=matrix(chpts,nrow=1)
    
#    if(msscan.debug.trace){
#         print("about to call Chisq.Contrib.C")
#         scan()
#    }

    chisq <- Chisq.Contrib.C(y, T, chpts) # chisq[i,j] holds contribution of sample j to segment i.

    yhat = matrix(0,nrow=N,ncol=T)
    sampleswithsegment=matrix(0,nrow=N,ncol=nrow(chpts))
    throw= rep(0,nrow(chpts))

    for(i in 1:nrow(chpts)){
        
        if(SAMPLE.CUTOFF.METHOD==1){
            sampleswithsegment[,i] = abs(apply(as.matrix(y[,chpts[i,1]:chpts[i,2]]),1,mean))>SAMPLE.HIGH.THRESH    
        }
        if(SAMPLE.CUTOFF.METHOD==2){
            sampleswithsegment[,i] = chisq[i,]>SAMPLE.THRESHOLD    
        }
        if(SAMPLE.CUTOFF.METHOD==3){
            sampleswithsegment[,i] = (abs(apply(as.matrix(y[,chpts[i,1]:chpts[i,2]]),1,mean))>SAMPLE.HIGH.THRESH) | (chisq[i,]>SAMPLE.THRESHOLD)
        }
        if(SAMPLE.CUTOFF.METHOD==4){
            sampleswithsegment[,i] = ((abs(apply(as.matrix(y[,chpts[i,1]:chpts[i,2]]),1,mean))>SAMPLE.LOW.THRESH) & (chisq[i,]>SAMPLE.THRESHOLD)) | (abs(apply(as.matrix(y[,chpts[i,1]:chpts[i,2]]),1,mean))>SAMPLE.HIGH.THRESH)
        }


        m2 <- y[which(sampleswithsegment[,i]==1),chpts[i,1]:chpts[i,2]]
        if( sum(sampleswithsegment[,i]==1) == 1 || length(chpts[i,1]:chpts[i,2]) == 1){
            m2 <- matrix(m2, nrow=sum(sampleswithsegment[,i]==1), ncol=length(chpts[i,1]:chpts[i,2])  )
        }
        yhat[which(sampleswithsegment[,i]==1),chpts[i,1]:chpts[i,2]] = apply(m2, 1, mean)
        if( sum(sampleswithsegment[,i])==0) throw[i] = 1
        
    }
    
    # Throw away those chpts i where sum(sampleswithsegment[,i])=0.
    chpts = chpts[which(!throw),]
    sampleswithsegment = sampleswithsegment[,which(!throw)]
    chisq = chisq[which(!throw),]
    
    list(yhat=yhat, chpts=chpts, sampleswithsegment=sampleswithsegment,chisq=chisq,Z=Z)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
