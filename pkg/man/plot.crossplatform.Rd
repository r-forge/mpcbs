\name{plot.crossplatform}
\alias{plot.crossplatform}

\title{Simultaneously plots DNA copy number data from multiple platforms. }
\description{
  This function plots DNA copy number data from multiple platforms in vertically aligned panels.  It aligns the position of probes across plots for easy comparison.  It can also
  show the positions of breakpoints in a joint segmentation.
}
\usage{
plot.crossplatform(pos, y, anchor, chpts = NULL, ranking = NULL, yhat = NULL, platform.names, num.panels = NA, xlim = NULL, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{pos}{A vector of sorted integer arrays, one array for each platform.  pos[[k]] should give the positions, in increasing order, of the probes of the k-th platform.}
  \item{y}{A vector of intensity levels for each platform.  y[[k]] should be the intensity levels corresponding to pos[[k]] for platform k.}
  \item{anchor}{The anchor set, returned by a call to merge.pos(...)}
  \item{chpts}{(optional) A set of positions of change-points, given in reference to the anchor set. }
  \item{ranking}{(optional) A set of integers showing the ranking of the change-points.}
  \item{yhat}{(optional)  A vector of estimated intensity levels, with yhat[[k]] being the fitted values for y[[k]].}
  \item{platform.names}{(optional) A list of platform names. }
  \item{num.panels}{(optional) The number of vertical planels to divide the plotting area (default=length(pos)).}
  \item{xlim}{ The limits of the x-axis in all of the plots. }
  \item{\dots}{ Other plotting arguments. }
}


\references{ Zhang, NR, Senbabaoglu, Y. and Li, J.Z. (2009) Joint Estimation of DNA Copy Number from Multiple Platforms.  Under review, download manuscript from
http://www-stat.stanford.edu/~nzhang/web_multiplatform/ }
\author{Nancy R. Zhang}
\seealso{ mpcbs, mpcbs.mbic, merge.pos }
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

# Plot the data

plot.crossplatform(pos,y,anchor=anchor,platform.names=platform.names, col="darkgray")

seg<-mpcbs.mbic(y,pos,anchor, MAX.CHPTS=50, platform.names=platform.names,plots=TRUE)
  
# Give it yhat:

plot.crossplatform(pos,y, yhat=seg$yhat, anchor=anchor, platform.names=platform.names, col="darkgray")

# Give it chpts, and it will show breaks:

plot.crossplatform(pos,y, yhat=seg$yhat, chpts=seg$chpts, anchor=anchor, platform.names=platform.names, col="darkgray")

# Very useful:  you can zoom in on the plots using parameter "xlim":

plot.crossplatform(pos,y, yhat=seg$yhat, chpts=seg$chpts, anchor=anchor, platform.names=platform.names, col="darkgray", xlim=c(5e7,1e8))

# Try zooming in further:

plot.crossplatform(pos,y, yhat=seg$yhat, chpts=seg$chpts, anchor=anchor, platform.names=platform.names, col="darkgray", xlim=c(5e7,6e7))



}


