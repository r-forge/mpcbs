\name{smoothing}
\alias{smoothing}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Smoothing}
\description{
This function uses Bounded Complexity MIXture (BCMIX) algorithm to estimate parent-specify copy numbers.
}
\usage{
smoothing(samplename, inputdata, platform = "Unknown", sm.plot = FALSE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{samplename}{a name indicates where the data is sampled from, eg. patient ID}
  \item{inputdata}{a file with four columns "Chromosome ID" (eg. any number from 1 to 23), "SNP Position (bp)", "logR" and "Bfrequency". The column names should be "Chr", "Position", "logR" and "bfreq". }
  \item{platform}{specifies where is data is from, eg. Illumina.}
  \item{sm.plot}{if TRUE, then the interim results will be ploted }
  \item{hyper}{a txt file specifies initial hyper parameters}
  \item{mixstate}{a txt file specifies initial states of SNPs}
  \item{niters}{!! not sure about its function }
  \item{nBIMIXiters}{specifies the number of loops for BIMIX }
  \item{K}{!! not sure about its function }
  \item{M}{!! not sure about its function }
  \item{outerloop}{specifies the number of loops for a complete run of BIMIX and state update }
  \item{selectHyper}{!! not sure about its function }
  \item{verbose}{if TRUE, print out process during the run }
  \item{save.interim}{if TRUE, interim cnvobj will be saved as .Rdata file}
}
\details{
    This function takes 2-3 minutes for a typical file containing logR and Bfrequency of SNP data from one human choromosome.
}
\value{
It returns files samplename.Chr1.Smooth.Rdata ect.
}
\references{ ~put references to the literature/web site here ~ }
\author{ ~~who you are~~ }
\note{ ~~further notes~~ 
}
\seealso{
this function follows with \code{\link{segmentation}} and \code{\link{pscnlist}}. 
\code{\link{pscn.plot}} is the plot function to view a variety of results.  
}
\examples{
data(ex0)  # this file contain four columns "Chromosome ID", "SNP Position (bp)", "logR" and "Bfrequency" with column names "Chr", "Position", "logR" and "bfreq". 
# smoothing(samplename="TCGA.02.0043.01A.01D", inputdata=mydata, platform="Illumina")  # to excute, remove the "#" at the begining
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
