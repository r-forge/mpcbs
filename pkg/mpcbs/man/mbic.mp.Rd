\name{mbic.mp}
\alias{mbic.mp}
\title{Computes the modified BIC for a given multi-platform segmentation. }
\description{
 Computes the modified BIC for a given multi-platform segmentation. 
}
\usage{
mbic.mp(S, imap, sigma, tau, rratio)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{S}{cumulative sums.}
  \item{imap}{imap object from an anchor set.}
  \item{sigma}{platform specific variances.}
  \item{tau}{a vector of ascending change-points. }
  \item{rratio}{response ratio for each platform. }
}
\value{
  \item{mbic}{mbic}
  \item{term1}{term 1 of the modified BIC}
  \item{term2}{term 2 of the modified BIC}
  \item{term3}{term 3 of the modified BIC}
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


}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{change-point}
\keyword{segmentation}
\keyword{DNA copy number}
\keyword{cross-platform integration}


