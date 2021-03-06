% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot_distributions.R
\name{plot_distributions}
\alias{plot_distributions}
\title{Plot the emppirical, binomial, and beta-binomial distributions}
\usage{
plot_distributions(minN, maxN, bins, empirical, d_combined_sorted_binned,
  e_combined_sorted_binned, yuplimit = 0.15)
}
\arguments{
\item{minN}{Integer. The minimum coverage level}

\item{maxN}{Integer. The maximum coverage level}

\item{bins}{Numeric. Breakpoints for the bins.}

\item{empirical}{The empirical distribution}

\item{d_combined_sorted_binned}{The weighted expected binomial distribution}

\item{e_combined_sorted_binned}{The weighted expected beta-binomial
distribution}
}
\description{
Plot the emppirical, binomial, and beta-binomial distributions
}
