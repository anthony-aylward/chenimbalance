#===============================================================================
# alleledb_calcOverdispersion.R
#===============================================================================

# Imports ======================================================================

#' @import VGAM




# Functions ====================================================================

#' @title Graded weights for SSE calculation
#'
#' @param r_min Range minimum
#' @param r_max Range maximum
#' @param bins Breakpoints mapped to unit interval
#' @return Graded weights for SSE calculation
#' @export
graded_weights_for_sse_calculation <- function(r_min, r_max, bins) {
  r <- seq(r_min, r_max, 2 * (r_max - r_min) / (length(bins) - 1))
  r <- r[2:length(r)]
  if ((length(bins)-1)%%2 != 0) {
    c(r, sort(r[1:(length(r)-1)], decreasing = TRUE)) 
  } else {
    c(r, sort(r[1:length(r)], decreasing = TRUE))
  }
}

#' @title Empirical allelic Ratio
#'
#' @details
#' Plot right closed interval left open (range] for x axis
#' note that you have pseudozeroes as counts in your data so thats fine
#'
#' @param data A data frame with columns "total" and "allelicRatio"
#' @param bins Breakpoints for bins
#' @return Histogram of allelic ratios
#' @export
empirical_allelic_ratio <- function(data, bins, maxN, minN = 6, plot = FALSE) {
  data.match <- data[data[["total"]] <= maxN & data[["total"]] >= minN,]
  h <- hist(
    data.match[["allelicRatio"]],
    xlim = range(0, 1),
    breaks = bins,
    right = TRUE,
    plot = plot
  )
  h[["counts"]] / sum(h[["counts"]])
}

#' Weight by empirical counts
#'
#' @param total Vector containing coverage per variant
#' @return A one-column matrix of weights for each coverage level
#' @seealso \code{\link{weighted_expected_binomial}}
#' @export
weight_by_empirical_counts <- function(total) {
  t <- as.data.frame(table(total), stringsAsFactors = FALSE)
  w <- matrix(0, max(total), 1)

  for (j in 1:nrow(t)) {
    w[as.integer(t[j, 1]), 1] <- t[j, 2]
  }
  w
}

#' Change "counts" into density
#'
#' @param d.combined.sorted.binned The combined, sorted, and binned
#'   pseudodensity
#' @return The combined, sorted, and binned density
#' @export
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
#' @param maxN Maximum coverage
#' @param p Binomial parameter (null probability)
#' @param w Data frame or matrix giving 
#' @return The weighted null beta/binomial probability mass function
#' @export
#' @seealso \code{\link{weighted_expected_binomial}}
#' @export
nulldistrib <- function(
  w,
  minN = 6,
  p = 0.5,
  binSize = 40,
  yuplimit = 0.15,
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
#' dbinom(seq(0,500),n=500,p=0.5) gives the pdf of 0-500, at n=500, p=0.5 in binomial distrib
#'
#' weight each probability by the number of SNPs with n reads
#'
#' @param total Vector containing coverage per variant
#' @param binSize Approximate number of bins
#' @return weighted expected binomial values
#' @export
weighted_expected_binomial <- function(
  total,
  minN = 6,
  p = 0.5,
  binSize = 40,
  yuplimit = 0.15,
  distrib = "binomial",
  b = 0
) {
  w <- weight_by_empirical_counts(total)
  d.combined.sorted.binned <- nulldistrib(
    w,
    minN = minN,
    p = p,
    binSize = binSize,
    yuplimit = yuplimit,
    distrib = distrib,
    b = b
  )
}