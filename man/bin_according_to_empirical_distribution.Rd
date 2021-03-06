% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/nulldistrib.R
\name{bin_according_to_empirical_distribution}
\alias{bin_according_to_empirical_distribution}
\title{Bin the combined and sorted pseudodistribution}
\usage{
bin_according_to_empirical_distribution(d.combined.sorted, binSize = 40)
}
\arguments{
\item{d.combined.sorted}{The combined and sorted pseudodistribution}

\item{binSize}{Integer. The approximate number of bins.}
}
\value{
The binned density
}
\description{
Bin it according to the empirical distribution
}
\details{
empirical right closed, left open ?hist; right=TRUE
(range] so no double counts
but zero gets excluded!!
}
\seealso{
\code{\link{nulldistrib}}
}
