\name{mpcbs}
\alias{mpcbs}
\title{Multi-platform circular binary segmentation with user tuned stoping threshold. }
\description{
    This function performs multi-platform circular binary segmentation for detecting DNA copy number changes.  The number of change-points is determined by a user given threshold.
}
\usage{
mpcbs(y, pos, anchors)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{y}{A vector of intensity levels for each platform.  y[[k]] should be the intensity levels corresponding to pos[[k]] for platform k.}
 \item{pos}{A vector of sorted integer arrays, one array for each platform.  pos[[k]] should give the positions, in increasing order, of the probes of the k-th platform.}
  \item{anchor}{The anchor set, returned by a call to merge.pos(...)}
   \item{MIN.SNPs}{The minimum number of SNPs for each platform in a segment. }
  \item{MIN.BP.LEN}{The minimum base pair length of a segment. }
  \item{WCHISQ.CUTOFF}{Cut-off for the projected chi-square statistic that determines when to stop the recursion. }
    \item{rratio}{ An array containing the signal response ratio of each platform.  }
  \item{platform.names}{ The names of the platforms. }
  \item{plots}{ A logical indicating whether to make progress plots. }
  \item{plotspdf}{ An optional pdf file where the progress plots will be recorded. }
  \item{use.filtered.scan}{ A logical indicating whether to use the faster filtered scan statistic (highly recommended). }

  
}
\details{
 MIN.SNPs is the minimum number of snps that a platform needs to have in a window to contribute to the combined scan statistic.  If the platform does not have enough SNPs, it will simply contribute 0 to the overall statistic, but if there is enough evidence from the other platforms, that window may still be called.
  
  MIN.BP.LEN is the minimum base pair length of any call. 
     
  If the rratio is not specified, it is assumed to be 1 for all platforms. 
  
  If use.filtered.scan is TRUE, the scan will roughly take O[nlog(n)] time.  If this is false, the scan will take O[n^2] time.
}
\value{

yhat=yhat,chpts=chpts.final,ranking=ranking,Z.score=maxZ

  \item{yhat }{the segmented (piece-wise constant) values for y.}
  \item{chpts }{the list of change-points.}
  
  \item{Z.score }{the z-score associated with each interval.}
  
  \item{ranking}{ranking of the intervals in terms of z-score.}
  ...
}\references{ Zhang, NR, Senbabaoglu, Y. and Li, J.Z. (2009) Joint Estimation of DNA Copy Number from Multiple Platforms.  Under review, download manuscript from
http://www-stat.stanford.edu/~nzhang/web_multiplatform/ }
\author{Nancy R. Zhang}
\seealso{ mpcbs, plot.crossplatform, merge.pos }

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
seg<-mpcbs.mbic(y,pos,anchor,WCHISQ.CUTOFF=10, platform.names=platform.names,plots=TRUE)


# Plot the data and segmentation, giving it yhat and change-points.
plot.crossplatform(pos,y, yhat=seg$yhat, chpts=seg$chpts, anchor=anchor, platform.names=platform.names, col="darkgray")


}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{change-point}
\keyword{segmentation}
\keyword{DNA copy number}
\keyword{cross-platform integration}


