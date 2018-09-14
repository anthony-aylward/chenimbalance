#===============================================================================
# optimize.R
#===============================================================================

#' @title optimize beta-binomial parameters
#'
#' @param total integer; total number of reads at each site
#' @param allelic_ratio numeric; allelic ratio at each site
#' @param minN minimum coverage level
#' @param binSize approximate number of bins
#' @return list
#' @export
alleledb_beta_binomial <- function(
  total,
  allelic_ratio,
  minN = 6,
  binSize = 40
) {
  bins <- pretty(0:1, binSize)
  maxN = min(2500, max(total))
  counts <- (
    data.frame(total = total, allelicRatio = allelic_ratio)[total >= minN,]
  )
  empirical <- empirical_allelic_ratio(
    counts,
    bins,
    maxN = maxN,
    minN = minN,
    plot = FALSE
  )
  w <- weight_by_empirical_counts(counts[["total"]])
  d_combined_sorted_binned <- nulldistrib(
    w,
    minN = minN,
    binSize = binSize
  )
  sse <- sum((empirical - d_combined_sorted_binned[,2])^2)
  w_grad <- graded_weights_for_sse_calculation(
    r_min = 0,
    r_max = 1,
    bins = bins
  )
  overdispersion_details <- choose_overdispersion_parameter(
    w_grad,
    w,
    empirical,
    sse
  )
  optimized_overdispersion_details <- optimize_overdispersion_parameter(
    w_grad,
    w,
    overdispersion_details[["b_and_sse"]],
    overdispersion_details[["b_choice"]],
    overdispersion_details[["sse"]],
    empirical,
    overdispersion_details[["counter"]],
    minN = minN,
    binSize = binSize
  )
  b = optimized_overdispersion_details[["b_choice"]]
  sse = optimized_overdispersion_details[["sse"]]
  prob_details <- choose_probability_of_success_parameter(
    w_grad,
    w,
    empirical,
    sse,
    b,
    r_by = 0.025
  )
  optimized_prob_details <- optimize_probability_of_success_parameter(
    w_grad,
    w,
    prob_details[["prob_and_sse"]],
    prob_details[["prob_choice"]],
    prob_details[["sse"]],
    empirical,
    prob_details[["counter"]],
    b,
    minN = minN,
    binSize = binSize
  )
  prob = optimized_prob_details[["prob_choice"]]
  sse = optimized_prob_details[["sse"]]
  list(
    prob = prob,
    b = b,
    sse = sse,
    shape1 = prob * (1 - b) / b,
    shape2 = (1 - prob) * (1 - b) / b
  )
}
