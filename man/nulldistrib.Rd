% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/nulldistrib.R
\name{nulldistrib}
\alias{nulldistrib}
\title{Weighted beta/binomial distribution}
\usage{
nulldistrib(w, minN = 6, p = 0.5, binSize = 40,
  distrib = "binomial", b = 0)
}
\arguments{
\item{w}{Data frame or matrix giving}

\item{minN}{Minimum coverage (min total num of reads (since it's left open right closed))}

\item{p}{Binomial parameter (null probability)}
}
\value{
The weighted null beta/binomial probability mass function
}
\description{
Weighted beta/binomial distribution
}
\details{
\code{d.combined} collects all results and correspond it to an allelic ratio:
col1=allelicRatio (based on binomial n=6, ar=0,1/6,2/6...)
col2=corresponding weighted value in binomial distribution, i.e.
pdf(n, k, p) * (num of empirical SNPs at n counts)
}
