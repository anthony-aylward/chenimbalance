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

#' @title bisect
#'
#' @param p p
#' @param p_sim p_sim
#' @param p_choice p_choice
#' @param fdr fdr
#' @param fdr_threshold fdr_threshold
#' @param by by
#' @param distrib distribution (binomial or betabinomial)
#' @param b overdispersion parameter
#' @param w w
#' @param p_thresh p_thresh
#' @return something
#' @export
bisect <- function(
  p,
  p_sim,
  p_choice,
  fdr,
  fdr_threshold,
  by,
  distrib = "binomial",
  b = 0,
  w,
  p_thresh
) {
  p_fdr_e <- matrix(0, 100, 3)
  e_prev <- 10
  flag <- 3
  ctr <- 1
  p_fdr_e[ctr, 1] <- p_choice
  p_fdr_e[ctr, 2] <- fdr
  p_fdr_e[ctr, 3] <- e_prev
  
  while (flag) {
    start <- max(0, (p_choice - by / 2))   
    end <- p_choice + by / 2
    by <- by / 4
    
    if (start == 0) {
      start <- 5e-4
    }
    
    range <- seq(start, end, by)
    
    for (i in range) {
      tp <- cutoff(i, p)
      
      if (distrib == "binomial") {
        fp <- fp(w,p_thresh, i, "binomial")
      } else if(distrib == "betabinomial") {
        fp <- fp(w, p_thresh, i, "betabinomial", b)
      }
      
      fdr_ind <- fp / tp
      e_curr <- fdr_threshold - fdr_ind
      ctr <- ctr + 1
      
      p_fdr_e[ctr, 1] <- i
      p_fdr_e[ctr, 2] <- fdr_ind
      p_fdr_e[ctr, 3] <- e_curr
      e_prev <- p_fdr_e[(ctr - 1), 3]
      p_choice <- i
      
      if (e_curr < 0){ break }
    }
    if (signif(p_fdr_e[ctr - 1, 3], 3) == signif(p_fdr_e[ctr, 3], 3)) {
      flag <- 0
    }
  }  
  return(p_fdr_e)
}
