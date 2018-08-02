#===============================================================================
# imbalance.R
#===============================================================================

#' @title FP
#'
#' @param w w
#' @param p Numeric. The binomial parameter.
#' @param p_thresh p_thresh
#' @return something
#' @export
fp <- function(
  w,
  p,
  p_thresh,
  distrib = "binomial",
  b = 0
) {
  a = lapply(as.integer(w[,1]), function(x) seq(0, x))
  
  if (distrib == "binomial") { 
    b = lapply(
        a,
        function(x) {
          apply(
            as.data.frame(2 * pbinom(x, max(x), p)),
            1,
            function(x) min(x, 1)
          )
        }
      )
  } else if (distrib == "betabinomial") { 
    b = lapply(
      a,
      function(x) {
        apply(
          as.data.frame(2 * pbetabinom(x, max(x), p, b)),
          1,
          function(x) min(x, 1)
        )
      }
    )
  }

  sum(
    sapply(
      mapply(
        function(x, y, z) x * y * z, 
        lapply(b, function(x) x <= p_thresh),
        b,
        w[,2]
      ),
      max
    )
  )
}

#' @title cutoff
#'
#' @param x The cutoff
#' @param y The argument
#' @return the number of values passing the cutoff
#' @export
cutoff <- function(x, y) {
  sum(y <= x)
}