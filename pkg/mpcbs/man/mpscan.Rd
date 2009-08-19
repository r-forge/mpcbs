\name{mpscan}
\alias{mpscan}
\title{Multi-platform scan using projected chi-square }
\description{
  Scans multi-platform data for all variable length windows have projected chi-square statistc higher than a fixed user-set threshold.
}
\usage{
mpscan(y, pos, anchors, win)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{y}{A vector of intensity levels for each platform.  y[[k]] should be the intensity levels corresponding to pos[[k]] for platform k.}
 \item{pos}{A vector of sorted integer arrays, one array for each platform.  pos[[k]] should give the positions, in increasing order, of the probes of the k-th platform.}
  \item{anchor}{The anchor set, returned by a call to merge.pos(...)}
  \item{win}{Largest allowed window size. }
  \item{MIN.SNPs}{Minimum number of SNPs allowed in a region.}
   \item{WCHISQ.CUTOFF}{Threshold level for the weighted chi-square statistic.}
}
\value{
 \item{yhat}{A vector of fitted y values, yhat[[k]] contains the fitted values to y[[k]].}
  \item{chpts}{A list of change-points, there are shared across platforms.}
 \item{platform.chisq}{The contribution of each platform to the overall chi-square value for each interval.}
 \item{Z}{The Z matrix.}
 
 }
\references{ Zhang, NR, Senbabaoglu, Y. and Li, J.Z. (2009) Joint Estimation of DNA Copy Number from Multiple Platforms.  Under review, download manuscript from
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


}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{change-point}
\keyword{segmentation}
\keyword{DNA copy number}
\keyword{cross-platform integration}


