#===============================================================================
# empirical_allelic_ratio.R
#===============================================================================

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
  print(names(data))
  data.match <- data[data[["total"]] <= maxN & data[["total"]] >= minN,]
  h <- hist(
    data.match[["allelicRatio"]],
    breaks = bins,
    right = TRUE,
    plot = plot
  )
  h[["counts"]] / sum(h[["counts"]])
}
