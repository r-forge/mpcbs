\name{cross.platform.consensus}
\alias{cross.platform.consensus}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Compute the cross platform consensus copy number given a segmentation. }
\description{
Estimate cross-platform consensus copy number given multi-platform copy number data and a segmentation.  This is done by minimizing a weighted least squares criterion. 
}
\usage{
cross.platform.consensus(segs, sigma2, platform.names, plots = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{segs}{ A listed of cross-platform segmentations, returned by mpcbs.mbic. }
  \item{sigma2}{ A computed variance for each platform (see example below).   }
  \item{platform.names}{ A list of names for each platform. }
  \item{plots}{ A logical indicating whether to make diagnostic plots. }
}

\value{
     \item{converged}{A logical indicating whether the weighted least squares has converged.}
     \item{consensus.cn}{A vector of same length as segs, consensus.cn[[k]] is the consensus copy number estimate over the anchor positions for the segmentation contained in segs[[k]].}
     \item{alpha}{The fitted platform specific offset.}
    \item{r}{The fitted platform specific response ratio.}
    \item{theta}{Fitted intensity for each platform in each segment.}
}

\references{ Zhang, NR, Senbabaoglu, Y. and Li, J.Z. (2009) Joint Estimation of DNA Copy Number from Multiple Platforms.  Under review, download manuscript from
http://www-stat.stanford.edu/~nzhang/web_multiplatform/ }
\author{Nancy R. Zhang}
\seealso{ mpcbs, mpcbs.mbic, plot.crossplatform }
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
seg<-mpcbs.mbic(y,pos,anchor, MAX.CHPTS=10, platform.names=platform.names,plots=TRUE)

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
xmin=5.4e7; xmax=5.7e7
plot.crossplatform(pos,y, yhat=seg$yhat, anchor=anchor,platform.names=platform.names, num.panels=K+1, xlim=c(xmin,xmax), col="darkgray")
plot(anchor$merged.pos, consensus.cn, type="l", lwd=2, col="red", 
        xlab="Genome Position", main="Consensus Copy Number Estimate", cex.lab=2, xlim=c(xmin,xmax))



}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{weighted least squares }

