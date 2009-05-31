\name{update.bic.baseline}
\alias{update.bic.baseline}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ ~~function to do ... ~~ }
\description{
  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
update.bic.baseline(this.y, segbag, carriers, this.sum.deltasq, this.sum.deltasqwin, newchpt, samp.ord, win, y.var, PRIOR.CARRIER.PROB = 0.1, PRIOR.SEGCOUNT = 10, do.plots = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{this.y}{ ~~Describe \code{this.y} here~~ }
  \item{segbag}{ ~~Describe \code{segbag} here~~ }
  \item{carriers}{ ~~Describe \code{carriers} here~~ }
  \item{this.sum.deltasq}{ ~~Describe \code{this.sum.deltasq} here~~ }
  \item{this.sum.deltasqwin}{ ~~Describe \code{this.sum.deltasqwin} here~~ }
  \item{newchpt}{ ~~Describe \code{newchpt} here~~ }
  \item{samp.ord}{ ~~Describe \code{samp.ord} here~~ }
  \item{win}{ ~~Describe \code{win} here~~ }
  \item{y.var}{ ~~Describe \code{y.var} here~~ }
  \item{PRIOR.CARRIER.PROB}{ ~~Describe \code{PRIOR.CARRIER.PROB} here~~ }
  \item{PRIOR.SEGCOUNT}{ ~~Describe \code{PRIOR.SEGCOUNT} here~~ }
  \item{do.plots}{ ~~Describe \code{do.plots} here~~ }
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
function(this.y,segbag,carriers,this.sum.deltasq,this.sum.deltasqwin,newchpt,samp.ord,win,y.var, 
                                PRIOR.CARRIER.PROB=0.1, PRIOR.SEGCOUNT = 10, do.plots=TRUE){
#    # for debugging:
#    this.y=t(y)          ##### remember to back transpose later.
#    newchpt=c(st,ed)
#    this.sum.deltasqwin=sum.deltasqwin[1:nrow(segbag)]
#    this.sum.deltasq = sum.deltasq[1:nrow(segbag)]
    
    
    T = nrow(this.y)
    N = ncol(this.y)
    
    m = nrow(segbag)+1    
    w = newchpt[2]-newchpt[1]+1
    
    # Each of these are length N vectors.
    if(newchpt[2]-newchpt[1]+1 >= 2){
        delta.hat = (apply(this.y[newchpt[1]:newchpt[2],],2,sum) - (w/T)*apply(this.y,2,sum))/(w*(1-w/T))
    } else {
     delta.hat = (this.y[newchpt[1]:newchpt[2],] - (w/T)*apply(this.y,2,sum))/(w*(1-w/T))
    }
    
    effect.hat = delta.hat^2/y.var
    sum.deltasq.loc = cumsum(effect.hat[samp.ord])  
    sum.deltasqwin.loc = cumsum(effect.hat[samp.ord]*w*(1-w/T))
    
    PRIOR.Nm = PRIOR.SEGCOUNT*N
    PRIOR.Carriers = PRIOR.CARRIER.PROB*PRIOR.Nm
    
    # Computing the prior probability of J, assuming a single phat=P(carrier) over all segments.
    sumJi = sum(carriers) + c(1:N) 
    Jprior.one.phat = sumJi*log((sumJi+PRIOR.Carriers)/(N*m+PRIOR.Nm)) + (N*m-sumJi)*log((N*m+PRIOR.Nm-sumJi-PRIOR.Carriers)/(N*m+PRIOR.Nm))
    
    # Computing the prior probability of J, assuming a different phat=P(carrier) over all segments.
    J = apply(carriers,2,sum)
    Jprior.sep.phat= sum( J*log((J+PRIOR.Carriers)/(N+PRIOR.Nm)) + (N-J)*log((N+PRIOR.Nm-J-PRIOR.Carriers)/(N+PRIOR.Nm)))
    Jprior.sep.phat = Jprior.sep.phat + c(1:N)*log((c(1:N)+PRIOR.Carriers)/(N+PRIOR.Nm)) + (N-c(1:N))*log((N+PRIOR.Nm-c(1:N)-PRIOR.Carriers)/(N+PRIOR.Nm))
         
    tauprior = - (lchoose(T,m) + log(m))    

    maxlik = 0.5*(sum(this.sum.deltasqwin) + sum.deltasqwin.loc)
    
    # do not take log if sum.deltasq, sum.deltasqwin is ZERO.
    
    betaprior.terms =  -sumJi/2 -(sumJi/2)*log((sum(this.sum.deltasqwin)+sum.deltasqwin.loc)/sumJi) 
    
    randomwalk.terms = 2*m*(1.6-1.5) - sum(log(this.sum.deltasq)) - log(sum.deltasq.loc)
    

    mbic.one.phat = maxlik + tauprior + Jprior.one.phat + betaprior.terms + randomwalk.terms 
    mbic.sep.phat = maxlik + tauprior + Jprior.sep.phat + betaprior.terms + randomwalk.terms 

    
    par(mfrow=c(3,1))
    plot(sqrt(effect.hat[samp.ord]*w*(1-w/T)),ylab="Chi")
    plot(maxlik,type="b", ylim=c(min(mbic.sep.phat, na.rm=TRUE),max(maxlik)))
    points(mbic.one.phat,type="b",col="red")
    segments(which.max(mbic.one.phat),0,which.max(mbic.one.phat), max(maxlik),col="red",lwd=1)
    points(mbic.sep.phat,type="b",col="blue")
    segments(which.max(mbic.sep.phat),0,which.max(mbic.sep.phat), max(maxlik),col="blue",lwd=1)
    plot(Jprior.one.phat,type="b",col="red", ylim=c(min(c(Jprior.sep.phat,betaprior.terms),na.rm=TRUE),0), ylab="Penalties",
        main=paste("PRIOR.SEGCOUNT=",PRIOR.SEGCOUNT,", PRIOR.CARRIER.PROB=",PRIOR.CARRIER.PROB))
    points(Jprior.sep.phat,type="b",col="blue")
    points(betaprior.terms,type="b",col="darkgreen")
    points(randomwalk.terms,type="b",col="orange")
     
    list(mbic.sep.phat = mbic.sep.phat, mbic.one.phat = mbic.one.phat, sum.deltasq.loc = sum.deltasq.loc, sum.deltasqwin.loc = sum.deltasqwin.loc)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
