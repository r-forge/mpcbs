\name{merge.pos}
\alias{merge.pos}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Merge the Probe Positions from Multiple Microarray Platforms.}
\description{
  Merges the probe positions across platforms, and gives a set of anchor points at which to test for change-points.  
  In other words, the anchor points are locations where the scan statistic for MPCBS is evaluated.  
}
\usage{
    merge.pos(pos, anchors = NULL)
    merge.pos(pos, anchors=1)
    
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{pos}{A vector of sorted integer arrays, one array for each platform.  pos[[k]] should give the positions, in increasing order, of the probes of the k-th platform.}
  \item{anchors}{A list of integers specifying which platform(s) to use as anchor.  If this is null, then the anchor is the union of the probe locations from all platforms. }
}
\details{
  It is necessary to call this function to obtain an anchor set before using any of the other functions in the mpcbs package.  In the scan max Z(s,t), the anchor points are the set of possibilities for (s,t).  Note that ALL probes from ALL platforms are used in computing Z(s,t), regardless of the set of anchor points used.  
  Smaller anchor set means less computing time. The default is to set the anchor to be the union of all probes on all platforms.    
}
\value{
  ~Describe the value returned
  If it is a LIST, use
  \item{merged.pos }{An increasing integer array containing the set of merged probe positions.}
  \item{imap }{Let N be the length of merged.pos, and K be the number of platforms.  This is an N by K matrix with imap(i,j) being the
  index of the first probe in the k-th platform that is strictly to the right of merged.pos[k], i.e. the index within pos[[k]] of the first value that is strictly larger than merged.pos[i].   }
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

anchor = merge.pos(pos)

}

