\name{compute.var}
\alias{compute.var}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Computes the variance of a piecewise constant time series that may contain jumps. }
\description{
Computes the variance of a piecewise constant time series that may contain jumps by using one-half the variance of the
lag-1 differences.  
}
\usage{
compute.var(y)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{y}{ A matrix of real numbers, each column is a piece-wise constant time series. }
  \item{use.mean}{If TRUE (default) the mean of the lag-1 variances is used, otherwise use the median. }
}

\value{
    a vector of variances, one for each column.
}
\author{Nancy R. Zhang}
\examples{
    a=c(rnorm(100), rnorm(100,1,1), rnorm(100,2,1))
    plot(a)
    compute.var(a)
}

