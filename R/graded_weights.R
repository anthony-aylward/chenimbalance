#===============================================================================
# graded_weights.R
#===============================================================================

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