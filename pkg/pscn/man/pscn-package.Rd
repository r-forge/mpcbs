\name{pscn-package}
\alias{pscn-package}
\alias{pscn}
\docType{package}
\title{
Parent Specific DNA Copy Number Inference
}
\description{
This package estimate parent specific DNA copy number for human genomes.
}
\details{
\tabular{ll}{
Package: \tab pscn\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2009-05-03\cr
License: \tab What license is it under?\cr
LazyLoad: \tab yes\cr
}
First prepare a file four columns "Chromosome ID" (eg. any number from 1 to 23), "SNP Position (bp)", "logR" and "Bfrequency". The column names should be "Chr", "Position", "logR" and "bfreq". Then use function \code{\link{smoothing}} to estimate parent-specify copy numbers. After that, use \code{\link{segmentation}} to determine change points and \code{\link{pscnlist}} to get a table showing major and minor copy numbers. \code{\link{pscn.plot}} is the plot function to view a variety of results.   } 
}
\author{
Who wrote it

Maintainer: Who to complain to <yourfault@somewhere.net>
~~ The author and/or maintainer of the package ~~
}
\references{
~~ Literature or other references for background information ~~
}
~~ Optionally other standard keywords, one per line, from file KEYWORDS in ~~
~~ the R documentation directory ~~
\keyword{ package }
\seealso{
\code{\link{smoothing}}, \code{\link{segmentation}}, \code{\link{pscnlist}}, \code{\link{pscn.plot}}
}
\examples{
# smoothing (estimate parent-specify copy numbers)
  data(ex0) # this file contain four columns "Chromosome ID", "SNP Position (bp)", "logR" and "Bfrequency" with column names "Chr", "Position", "logR" and "bfreq".
  smoothing(samplename="TCGA.02.0043.01A.01D", inputdata=mydata, platform="Illumina") 

# segmentation (determine change points) this function can only be run after the function "smoothing".
  segmentation(samplename="TCGA.02.0043.01A.01D",chrs=1,regroup.percent=0.05) 

# pscnlist (get a table showing major and minor copy numbers) this function can only be run after the function "segmentation".
  pscnlist(samplename="TCGA.02.0043.01A.01D", chrs=1)

# pscn.plot (more example at \link{pscn.plot}) Except the plot for "copy number", which can only be done after the function "pscnlist", other plots can be done as soon as the function "segmentation" is finished.
  mylist = read.table("TCGA.02.0043.01A.01D.shortlist.txt",header=TRUE) # read in a Text document named samplename.shortlist.txt for "copy number" plot 
  pscn.plot(cndata=mylist, which.plot="copy number", chrid=1)
  load("TCGA.02.0043.01A.01D.Chr1.pscn.Rdata") # load samplename.Chr1.pscn.Rdata  (Chr1 can be replaced by any other chromosomes, such as Chr2, etc.) for plot other than "copy number".
  pscn.plot(cnvobj, cndata=mylist, which.plot="bfreq", changepoint=FALSE)
}
