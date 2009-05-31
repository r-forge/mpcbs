\name{fscan.max}
\alias{fscan.max}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ ~~function to do ... ~~ }
\description{
  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
fscan.max(this.S, this.SST, y.var = NULL, use.BY.statistic = FALSE, use.Project.statistic = FALSE, this.imap = NULL, f = NULL, MIN.SNPs = 2, ALPHA = 0, delta = c(1, -1), rratio = NULL, doplots = FALSE, verbose = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{this.S}{ ~~Describe \code{this.S} here~~ }
  \item{this.SST}{ ~~Describe \code{this.SST} here~~ }
  \item{y.var}{ ~~Describe \code{y.var} here~~ }
  \item{use.BY.statistic}{ ~~Describe \code{use.BY.statistic} here~~ }
  \item{use.Project.statistic}{ ~~Describe \code{use.Project.statistic} here~~ }
  \item{this.imap}{ ~~Describe \code{this.imap} here~~ }
  \item{f}{ ~~Describe \code{f} here~~ }
  \item{MIN.SNPs}{ ~~Describe \code{MIN.SNPs} here~~ }
  \item{ALPHA}{ ~~Describe \code{ALPHA} here~~ }
  \item{delta}{ ~~Describe \code{delta} here~~ }
  \item{rratio}{ ~~Describe \code{rratio} here~~ }
  \item{doplots}{ ~~Describe \code{doplots} here~~ }
  \item{verbose}{ ~~Describe \code{verbose} here~~ }
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
function(this.S, this.SST,  y.var=NULL, use.BY.statistic=FALSE, use.Project.statistic=FALSE, this.imap=NULL, 
                    f=NULL,MIN.SNPs=2,ALPHA=0,delta=c(1,-1), rratio=NULL, doplots=FALSE, verbose=FALSE){
  T=nrow(this.S)   # Number of SNPs
  N=ncol(this.S)   # Number of samples
  if(T<=MIN.SNPs){
    maxZ=NA
    seg=c(NA,NA)
  } else {

      if(is.null(f)){
        f.power = seq(min(-floor(log10(T))+1,0),0,1)
        f=10^f.power
      }
      
      if(T<1000) f=1 # even if a value for f is passed in, do not allow filtering for short sequences.
      
      f = sort(f)
      if(f[length(f)] != 1) f = c(f,1)
      R = length(f)
      L = ceiling(f[2:R]/f[1:(R-1)])
      L = c(ceiling(T*f[1]),L)
      
      chpts = matrix(ncol=2,nrow=0)
      chpts.Z = matrix(nrow=1,ncol=0)
    
      if(is.null(this.imap)){
        this.imap = matrix(rep(seq(1,T,1),N),ncol=N,byrow=FALSE)
      }
    
      g2<-function(u){ u^2*exp(u^2/2)/(ALPHA+exp(u^2/2))}
      g2.moments=computeMoments(g2,0);
      if(is.null(rratio)) rratio = rep(1,N)
      
      for(r in 1:R){
        stepsize = floor(1/f[r])
        
        t = seq(1,T,stepsize) # t is the filtered anchor set.
        if(t[length(t)]<T) t = c(t,T) # always include the last datapoint in the set.
    
        f.S = this.S[t,]
        f.imap = this.imap[t,]
    
        # produce a diagnostic plot.
        if(doplots){
          plot(this.S[,1],xlim=c(2000,3000),xlab="SNP #", ylab="Cumulative sum of sample 1", main=paste("Round ",r,": Take 1 out of every ", stepsize," points.",sep=""))
          lines(t,f.S[,1],col="red")
          points(t,f.S[,1],col="red",pch=17,cex=2)
        }
        
        # Refine previously found change-points using the denser anchor set.
        if(nrow(chpts)>0){
          for(i in 1:nrow(chpts)){
            ind.L = chpts[i,1] \%/\% stepsize
            ind.R = chpts[i,2] \%/\% stepsize
            
            check.win = f[r]/f[r-1]
            start.inds = c((ind.L - check.win):(ind.L+check.win))  
            end.inds = c((ind.R-check.win):(ind.R+check.win))
            
                  
            if(use.BY.statistic){
                Z.part=ComputeBYZ.fromS.R.partial(f.S,y.var,f.imap,start.inds, end.inds,delta,MIN.SNPS)
                Z.part = Z.part/N
            } else {
                if(use.Project.statistic){
                    Z.part=ComputeProjectedZ.fromS.R.partial(f.S,this.SST,f.imap,start.inds, end.inds,rratio,MIN.SNPs)
                    Z.part = (Z.part -1)/sqrt(2)          
                } else {
                    Z.part=ComputeZ.fromS.R.partial(f.S,this.SST,f.imap,start.inds, end.inds,ALPHA,MIN.SNPs)
                    Z.part = (Z.part - N*g2.moments$psidot)/sqrt(g2.moments$psidotdot*N)          
                }
            }
    
            #  image.plot(start.inds,end.inds,Z.part, xlab="Start of change - 1", ylab="End of change")
    
            maxind = matrix.max(Z.part)
            improved.cp = c(t[start.inds[maxind[1]]]+1, t[end.inds[maxind[2]]])
            improved.Z = Z.part[maxind[1],maxind[2]]
            if(verbose) cat("fscan.max: Changepoints (",chpts[i,1],", ",chpts[i,2],") refined to (", improved.cp[1],", ",improved.cp[2],").  Z-score improved from ", chpts.Z[i]," to ",improved.Z,"\n", sep="")
            chpts[i,] = improved.cp
            chpts.Z[i] = improved.Z
          }
        }
        
        
        # Do scan with min window size MIN.SNPS and max window size L_i,    
        if(use.BY.statistic){
            Z = ComputeBYZ.fromS.R(f.S, y.var, f.imap, L[r], delta, MIN.SNPS)
            Z = Z/N
        } else {
            if(use.Project.statistic){
                Z = ComputeProjectedZ.fromS.R(f.S, this.SST, f.imap, L[r], rratio, MIN.SNPs)
                Z = (Z - 1)/sqrt(2)    
            } else {
                # Do scan with min window size MIN.SNPs and max window size L_i,
                Z = ComputeZ.fromS.R(f.S, this.SST, f.imap, L[r], ALPHA, MIN.SNPs)
                Z = (Z - N*g2.moments$psidot)/sqrt(g2.moments$psidotdot*N)      
            }
        }
        # produce a diagnostic plot.
        if(doplots) image.plot(t,c(1:L[r]),Z, xlab="Start of change - 1", ylab="Window size", main=paste("Sum of chisquare (Z-score) for round ",r,sep=""))
    
        maxind = matrix.max(Z)
        #  newchpt= c(t[maxind[1]],t[maxind[1]+maxind[2]])
        newchpt= c(t[maxind[1]],t[maxind[1]+maxind[2]-1])  # changed 10/22, if maxind[1]+maxind[2] might become length(t)+1...
        newZ = Z[maxind[1],maxind[2]]
        if(verbose){
          cat("fscan.max: Change-point found in round ",r,":\n")
          print(c(newchpt,newZ), digits=1)
        }
        chpts = rbind(chpts,newchpt)
        if(length(chpts) == 2) chpts =matrix(nrow=1,ncol=2,data=chpts)      
        chpts.Z = c(chpts.Z,newZ)
      }
    
      # Take the maximum over the rounds, return chpt and Z score.
      max.ind = which.max(chpts.Z)
      maxZ = chpts.Z[max.ind]
      seg = chpts[max.ind,]
  }
  list(maxZ = maxZ, seg = seg)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
