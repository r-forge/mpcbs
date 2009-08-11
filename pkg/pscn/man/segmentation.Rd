\name{segmentation}
\alias{segmentation}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Segmentation}
\description{
This function determines change points after the smoothing function.
}
\usage{
segmentation(samplename, chrs = 1:23, gender = NULL, MIN.SNPS = 20, regroup.percent = -1, combine.alpha = 0.05, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{samplename}{ a name indicates where the data is sampled from, eg. patient ID }
  \item{chrs}{indicates which chromosomes of the sample are being segmented}
  \item{gender}{This function will determine the gender of the person where the data comes from if 23 is included in the "chrs" no match what input of "gender" is.  However, if 23 is not included in the chrs, the gender of the person where data comes from will be determined by the gender input.}
  \item{MIN.SNPS}{specifies the minimum number of SNPs in a segment.}
  \item{MIN.HET.SNPS}{specifies the minimum number of heterozyous SNPs in a segment.} 
  \item{regroup.percent}{when regroup using Bfrequency, if the obtained percent of heterozyous SNPs is less than the regroup.percent, some current homozyous SNPs will be classified as heterozyous.  The default value here is "-1", that is, this criterion is not used as default. If would like to use it, can specify some value, such as "0.2".}
  \item{combine.alpha}{the significance level of hypotheses tesing for combining adjacent segments. }
  \item{verbose}{if TRUE, print out process during the run}
  \item{DALPHA.THRESH}{The critical value for being a change point if the $\alpha$ value of the current SNP is
different from the previous SNP by DALPHA.THRESH.}
  \item{DBETA.THRESH}{The critical value for being a change point if the $\beta$ value of the current SNP is
different from the previous SNP by DBETA.THRESH.}
  \item{DELTATHRESH}{Delete a change point if its SNP index increases less than DELTATHRESH from the previous change point.}
  \item{ALPHA.THRESH, BETA.THRESH}{used when merge=TURE. Combine two neighbor segments when
the difference of alpha values between the two segments is smaller than ALPHA.THRESH and
the difference of beta values between the two segments is smaller than BETA.THRESH.}
  \item{BASELINE.BETA.THRESH}{set beta value to 0 if it is less than BASELINE.BETA.THRESH.}
  \item{FRAC.CHPTS}{used to calculate DALPHA.THRESH and BETA.THRESH if they are not specified.}
}
\details{
This function first uses striengent criteria to find change points. Then these roughly found segments are regrouped using Bfreqency because in each segment, there are at most four clusters of Bfreqence. After that, adjacent segments are tested to determine whether they can be combined.  
}
\value{
It returns files samplename.Chr1.Segment.Rdata, samplename.Chr1.pscn.Rdata, ect.
}
\references{ ~put references to the literature/web site here ~ }
\author{ ~~who you are~~ }
\note{follows the result of \link{smoothing}.  This function can only be run after the function \link{smoothing}.}
\seealso{
Raw logR and Bfreqency data should use be analyzed using \link{smoothing} first, then follows this function. \link{pscnlist} follow this function to generate allele specific copy numbers.  \link{pscn.plot} is the plot function to view a variety of results.  
}
\examples{
# follows the result of \link{smoothing}.  This function can only be run after the function \link{smoothing}.
segmentation(samplename="TCGA.02.0043.01A.01D",chrs=1) 
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
