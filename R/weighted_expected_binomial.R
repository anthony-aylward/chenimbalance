#===============================================================================
# weighted_expected_binomial.R
#===============================================================================

# Imports ======================================================================

#' @import VGAM




# Functions ====================================================================

#' @title Change "counts" into density
#'
#' @param d.combined.sorted.binned The combined, sorted, and binned
#'   pseudodensity
#' @return The combined, sorted, and binned density
#' @export
#' @seealso \code{\link{bin_according_to_empirical_distribution}}
change_counts_into_density <- function(d.combined.sorted.binned) {
  d.combined.sorted.binned[,2] <- d.combined.sorted.binned[,2] / sum(
    d.combined.sorted.binned[,2]
  )
  d.combined.sorted.binned
}

#' Bin the combined and sorted pseudodistribution
#'
#' Bin it according to the empirical distribution
#'
#' empirical right closed, left open ?hist; right=TRUE
#' (range] so no double counts
#' but zero gets excluded!!
#'
#' @param d.combined.sorted The combined and sorted pseudodistribution
#' @param binSize Integer. The approximate number of bins.
#' @return The binned density
#' @export
#' @seealso \code{\link{nulldistrib}}
bin_according_to_empirical_distribution <- function(
  d.combined.sorted,
  binSize = 40
) {
  bins <- pretty(0:1, binSize)
  start <- 0
  end <- 0
  d.combined.sorted.binned <- matrix(0, length(bins) - 1, 2)
  
  for (z in 2:length(bins)) { ##skip 0
    start <- bins[z - 1]
    end <- bins[z]
    row <- z - 1
    d.combined.sorted.binned[row, 1] <- mean(c(end, start))
    
    d.combined.sorted.binned[row, 2] <-  sum(
      d.combined.sorted[
        (d.combined.sorted[,1] <= end & d.combined.sorted[,1] > start),
        2
      ]
    )
    
    if (row == 1) {
      d.combined.sorted.binned[row, 2] = sum(
        d.combined.sorted[
          (d.combined.sorted[,1] <= end & d.combined.sorted[,1] >= start),
          2
        ]
      )
    }    
  }

  change_counts_into_density(d.combined.sorted.binned)
}

#' Weight each coverage level with actual counts in empirical
#'
#' @param d.combined The combined pseudodistribution
#' @param d The parametric density value
#' @param w w
#' @param k k
#' @param minN Integer. The minimum coverage level.
#' @return The combined pseudodistribution weighted by the empirical
#'   distribution
#' @export
#' @seealso \code{\link{nulldistrib}}
weight_with_actual_counts_in_empirical <- function(
  d.combined,
  d,
  w,
  i,
  k,
  ptr,
  minN = 6
) {
  d.w <- d * w[i, 1]
    
  if (i == minN) {
    d.combined[ptr:length(k), 1] <- k / i
    d.combined[ptr:length(k), 2] <- d.w
    colnames(d.combined) = c('allelicRatio', 'wBinDist')
  } else {
    d.combined[ptr:(ptr + length(k) - 1), 1] <- k / i
    d.combined[ptr:(ptr + length(k) - 1), 2] <- d.w
  }
  d.combined
}

#' @title A parametric probability mass value
#'
#' @description From a binomial or beta-binomial distribution
#'
#' @param k k
#' @param i i
#' @param distrib Character string. The distribution to draw from, either
#'   "binomial" or "betabinomial"
#' @param p Numeric. The binomial parameter.
#' @param b Numeric. The beta parameter.
#' @export
#' @seealso \code{\link{nulldistrib}}
parametric_probability_mass <- function(
  k,
  i,
  distrib = "binomial",
  p = 0.5,
  b = 0
) {
  if (distrib == "binomial") {
    dbinom(k, i, p)
  } else if (distrib == "betabinomial") {
    dbetabinom(k, i, p, b)
  }
}

#' @title Weighted beta/binomial distribution
#'
#' @details
#' \code{d.combined} collects all results and correspond it to an allelic ratio:
#' col1=allelicRatio (based on binomial n=6, ar=0,1/6,2/6...)
#' col2=corresponding weighted value in binomial distribution, i.e.
#' pdf(n, k, p) * (num of empirical SNPs at n counts)
#'
#' @param minN Minimum coverage (min total num of reads (since it's left open right closed))
#' @param p Binomial parameter (null probability)
#' @param w Data frame or matrix giving 
#' @return The weighted null beta/binomial probability mass function
#' @export
#' @seealso \code{\link{weighted_expected_binomial}}
nulldistrib <- function(
  w,
  minN = 6,
  p = 0.5,
  binSize = 40,
  distrib = "binomial",
  b = 0
) {
  maxN <- min(2500, nrow(w))
  d.combined <- matrix(0, sum(seq(minN + 1, maxN + 1)), 2)
  ptr <- 1
  for (i in minN:maxN) {
    k <- seq(0, i)
    d <- parametric_probability_mass(k, i, distrib = distrib, p = p, b = b)
    d.combined <- weight_with_actual_counts_in_empirical(
      d.combined,
      d,
      w,
      i,
      k,
      ptr,
      minN = minN
    )
    ptr <- ptr + length(k)
  }
  d.combined.sorted <- d.combined[order(d.combined[,1], d.combined[,2]),]
  bin_according_to_empirical_distribution(d.combined.sorted, binSize = binSize)
}

#' Weighted expected binomial 
#'
#' Plot the binomial distribution for each n
#'
#' dbinom(seq(0,500),n=500,p=0.5) gives the pdf of 0-500, at n=500, p=0.5 in
#' binomial distrib
#'
#' weight each probability by the number of SNPs with n reads
#'
#' @param total Vector containing coverage per variant
#' @param binSize Approximate number of bins
#' @return weighted expected binomial values
#' @export
weighted_expected_binomial <- function(
  w,
  minN = 6,
  p = 0.5,
  binSize = 40,
  distrib = "binomial",
  b = 0
) {
  w <- weight_by_empirical_counts(total)
  d.combined.sorted.binned <- nulldistrib(
    w,
    minN = minN,
    p = p,
    binSize = binSize,
    distrib = distrib,
    b = b
  )
}