#===============================================================================
# plot_distributions.R
#===============================================================================

#' @title Plot the emppirical, binomial, and beta-binomial distributions
#'
#' @param minN Integer. The minimum coverage level
#' @param maxN Integer. The maximum coverage level
#' @param bins Numeric. Breakpoints for the bins.
#' @param empirical The empirical distribution
#' @param d_combined_sorted_binned The weighted expected binomial distribution
#' @param e_combined_sorted_binned The weighted expected beta-binomial
#'   distribution
#' @export
plot_distributions <- function(
  minN,
  maxN,
  bins,
  empirical,
  d_combined_sorted_binned,
  e_combined_sorted_binned,
  yuplimit = 0.15
) {
  barplot(
    empirical,
    ylim = c(0, yuplimit),
    ylab = "density",
    xlab = "allelicRatio",
    names.arg = bins[2:length(bins)] - bins[[2]] / 2,
    main = paste("n=", minN, "-", maxN)
  )
  par(new = TRUE)
  plot(
    d_combined_sorted_binned,
    ylim = c(0, yuplimit),
    pch = 16,
    type = "b",
    col = "red",
    bty = "n",
    ylab = "",
    xlab = "",
    yaxt = "n",
    xaxt = "n",
    yaxs = "i"
  )
  par(new = TRUE)
  plot(
    e_combined_sorted_binned,
    ylim = c(0, yuplimit),
    pch = 16,
    type = "b",
    col = "blue",
    bty = "n",
    ylab = "",
    xlab = "",
    yaxt = "n",
    xaxt = "n",
    yaxs = "i"
  )
}