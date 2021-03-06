% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/overdispersion.R
\name{choose_overdispersion_parameter}
\alias{choose_overdispersion_parameter}
\title{Choose overdispersion parameter}
\usage{
choose_overdispersion_parameter(w_grad, w, empirical, sse, minN = 6,
  binSize = 40, b_choice = 0, r_sta = 0.1, r_end = 0.99,
  r_by = 0.1, p = 0.5, cores = detectCores())
}
\arguments{
\item{w_grad}{w_grad}

\item{w}{w}

\item{empirical}{empirical}

\item{sse}{sse}

\item{minN}{minimum coverage level}

\item{binSize}{approximate number of bins}

\item{r_sta, r_end}{start and end of the parameter range}

\item{r_by}{step of the parameter range}

\item{cores}{number of cores to use}
}
\value{
list
}
\description{
Minimize sse for betabinomials
}
