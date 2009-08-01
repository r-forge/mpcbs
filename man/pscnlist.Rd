\name{pscnlist}
\alias{pscnlist}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Generate a List of Allele Specific Copy Numbers}
\description{
This function generates a list of allele specific copy numbers of a sample
}
\usage{
pscnlist(samplename, chrs = 1:23, MIN.SNPS = 20, alpha = 0.01, FWalpha = 0.01, threshold = 0.18, hard.threshold=0.125, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{samplename}{a name indicates where the data is sampled from, eg. patient ID }
  \item{chrs}{indicates which chromosomes of the sample are being segmented}
  \item{MIN.SNPS}{specifies the minimum number of SNPs in a segment. }
  \item{alpha}{significance level for each segment}
  \item{FWalpha}{ significance level for the family of segments with changed copy numbers }
  \item{threshold}{In a segment, if copy numbers for both chromosomes are considered to be changed (hypothesis tests significant), then the one having copy number closer to the normal copy number will be tested to see if it differs from the normal copy number as large as this threshold; if not, treat it as normal.}
  \item{hard.threshold}{In a segment, if copy numbers for both chromosomes are considered to be changed (hypothesis tests significant).  Then, both will be tested to see if they differ from the normal copy number as large as this threshold; if both not, consider both to be normal.  This threshold should be smaller than the "threshold" above}
  \item{cutvalue}{stop the iteration of estimating major and minor copy numbers if the difference in L2-norm is less than this cut value}
  \item{pi1}{mixed Gaussian parameter, initial percent of one Gaussian}
  \item{mu1.initial}{mixed Gaussian parameter, initial mean of one Gaussian}
  \item{mu2.initial}{mixed Gaussian parameter, initial mean of the other Gaussian}
  \item{sigma1.initial}{mixed Gaussian parameter, initial variance of one Gaussian}
  \item{sigma2.initial}{mixed Gaussian parameter, initial variance of the other Gaussian}
  \item{total0}{For each segment, a test is carried out using a modified BIC criterion to determine the number of clusters the heterozyous SNPs. If the number of cluster is 1, then the total copy number is compared with this value to determine whether this segment has both chromosomes deleted.}
  \item{total1}{similarly usage as total0, while used for determining if the segment has both chromosomes increased if the total copy number is larger than total1} 
}
\details{
In this step, the segments with no change are first determined using a modified BIC criterion, then the normal mean and variance.  For other segments, a mixed Gaussian model is used to estimate the mean and variance of the two component Gaussian densities.  They are compared with the normal mean and variance to determine whether its copy number changes.
}
\value{
It returns three .txt files named samplename.longlist.txt, samplename.longlist.update.txt, samplename.shortlist.txt 
}
\references{ ~put references to the literature/web site here ~ }
\author{ ~~who you are~~ }
\note{follows the result of \link{segmentation}.  this function can only be run after the function \link{segmentation}.}
\seealso{
 Raw logR and Bfreqency data should use be analyzed using \link{smoothing} and \link{segmentation} first, then follows this function.  \link{pscn.plot} is the plot function to view a variety of results.  
}
\examples{
# follows the result of \link{segmentation}.  this function can only be run after the function \link{segmentation}.
  pscnlist(samplename="TCGA.02.0043.01A.01D", chrs=1)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
