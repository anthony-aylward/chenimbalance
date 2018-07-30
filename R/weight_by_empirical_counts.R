#===============================================================================
# weight_by_empirical_counts.R
#===============================================================================

#' @title Weight by empirical counts
#'
#' @param total Vector containing coverage per variant
#' @return A one-column matrix of weights for each coverage level
#' @seealso \code{\link{weighted_expected_binomial}}
#' @export
#' @seealso \code{link{weighted_expected_binomial}},
#'   \code{link{choose_overdispersion_parameter}}
weight_by_empirical_counts <- function(total) {
  t <- as.data.frame(table(total), stringsAsFactors = FALSE)
  w <- matrix(0, max(total), 1)
  for (j in 1:nrow(t)) {
    w[as.integer(t[j, 1]), 1] <- t[j, 2]
  }
  w
}