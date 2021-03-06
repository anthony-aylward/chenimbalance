% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/nulldistrib.R
\name{parametric_probability_mass}
\alias{parametric_probability_mass}
\title{A parametric probability mass value}
\usage{
parametric_probability_mass(k, i, distrib = "binomial", p = 0.5,
  b = 0)
}
\arguments{
\item{k}{k}

\item{i}{i}

\item{distrib}{Character string. The distribution to draw from, either
"binomial" or "betabinomial"}

\item{p}{Numeric. The binomial parameter.}

\item{b}{Numeric. The beta parameter.}
}
\description{
From a binomial or beta-binomial distribution
}
\seealso{
\code{\link{nulldistrib}}
}
