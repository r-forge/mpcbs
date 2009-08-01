\name{mpcbs-package}
\alias{mpcbs-package}
\docType{package}
\title{
Multi-platform circular binary segmentation for analysis of DNA copy number data.  
}
\description{
    This library contains functions that perform multi-platform circular binary segmentation for detecting DNA copy number changes using 
    intensity data from multiple platforms or technical replictes from the same platform.  
    The fineness of the segmnetation (i.e. the number of change-points) can be determined either using a modified
    BIC criterion or a pre-set weighted chi-square threshold.  
}
\details{
\tabular{ll}{
Package: \tab mpcbs\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2009-07-23\cr
License: \tab Academic\cr
LazyLoad: \tab yes\cr
}
The main segmentation functions are "mpcbs" and "mpcbs.mbic".   
The function "mpcbs" does multiplatform segmentation with user defined
stopping criterion, and the function "mpcbs.mbic" does multiplatform 
segmentation with stopping criterion determined by the BIC criterion.
After segmentation, the function "cross.platform.consensus" computes the platform 
specific response ratios and gives a cross-platform consensus copy number estimate
for each segment.  These two functions ("mpcbs"/"mpcbs.mbic" and "cross.platform.consensus")
can be iterated for better performance.  

Several plotting functions for visualizing multi-platform data are also provided.  The most useful
of these is plott.crossplatform.  Also, scatter.ci and scatter.cpfit are useful for diagnosing
the fit given by cross.platform.consensus.

An additional function, mpscan, does not do circular binary segmentation but simply does a sliding
window scan, with varying window width, using the projected chi-square statistic.

}
\author{
Nancy R. Zhang <nrzhang@stanford.edu>
}
\references{
Zhang, NR, Senbabaoglu, Y. and Li, J.Z. (2009) Joint Estimation of DNA Copy Number from Multiple Platforms.  Under review, download manuscript from
http://www-stat.stanford.edu/~nzhang/web_multiplatform/
}
\keyword{ package }
\keyword{change-point}
\keyword{segmentation}
\keyword{DNA copy number}
\keyword{cross-platform integration}

\seealso{
    msscan
}
\examples{

data(mpcbs.example)

# there are 3 platforms represented in this data example: Illumina, Affymetrix, and Agilent.
names(mpcbs.example)

# K is the number of platforms.
K=3

# Store the chromosome positions in vector pos,
# the intensities in vector y:

pos=vector("list",K)
pos[[1]] = mpcbs.example$illu[,1]
pos[[2]] = mpcbs.example$affy[,1]
pos[[3]] = mpcbs.example$agil[,1]

y = vector("list",K)
y[[1]] = mpcbs.example$illu[,2]
y[[2]] = mpcbs.example$affy[,2]
y[[3]] = mpcbs.example$agil[,2]

# Names of the platforms:
platform.names=c("Illumina","Affymetrix","Agilent")

# Get the anchor set.
anchor = merge.pos(pos)

# Perform the segmentation.
seg<-mpcbs.mbic(y,pos,anchor, MAX.CHPTS=50, platform.names=platform.names,plots=TRUE)

# Compute platform specific variance terms.
sigma2 = rep(0,K)
for(k in 1:K) sigma2[k] = compute.var(y[[k]])
# When data from multiple chromosomes are available,
# put them in the list "segs".  Now, the list has length 1.
segs=vector("list",1)
segs[[1]] = seg
res=cross.platform.consensus(segs, sigma2, platform.names, plots=TRUE)
consensus.cn = res$consensus.cn[[1]]    # consensus.cn[t] is consensus copy number estimate at anchor$merged.pos[t]

# Finally, plot the cross platform consensus.

# Set parameter num.panels to K+1 (default is K), because we will plot consensus in (K+1)th panel.
plot.crossplatform(pos,y, yhat=seg$yhat, anchor=anchor,platform.names=platform.names, num.panels=K+1, col="darkgray")
plot(anchor$merged.pos, consensus.cn, type="l", lwd=2, col="red", 
        xlab="Genome Position", ylab="Log Ratio", main="Consensus Copy Number Estimate", cex.lab=2)
    
# Zoom in:
xmin=5e7; xmax=1e8
plot.crossplatform(pos,y, yhat=seg$yhat, anchor=anchor,platform.names=platform.names, num.panels=K+1, xlim=c(xmin,xmax), col="darkgray")
plot(anchor$merged.pos, consensus.cn, type="l", lwd=2, col="red", 
        xlab="Genome Position", ylab="Log Ratio", main="Consensus Copy Number Estimate", cex.lab=2, xlim=c(xmin,xmax))


}
