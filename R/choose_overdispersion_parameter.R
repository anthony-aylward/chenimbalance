#===============================================================================
# choose_overdispersion_parameter.R
#===============================================================================

#' Choose overdispersion parameter
#'
#' Minimize sse for betabinomials
#'
#' @param data A data frame with columns "total" and "allelicRatio"
#' @param bins Breakpoints for bins
#' @return Histogram of allelic ratios
#' @export
choose_overdispersion_parameter <- function(
  w_grad,
  w,
  empirical,
  sse,
  minN = 6,
  binSize = 40,
  b_choice = 0,
  r_sta = 0.1,
  r_end = 0.99,
  r_by  = 0.1
) {
  counter <- 1
  b_and_sse = matrix(
    c(b_choice, sse, rep(0, 48)),
    nrow = 50,
    ncol = 2,
    byrow = TRUE,
    dimnames = list(NULL,  c("b", "sse"))
  )
  labels <- matrix(0, nrow = 50, ncol = 1)
  b_range <- seq(r_sta, r_end, by = r_by)
  for (k in b_range) {
    e_combined_sorted_binned <- nulldistrib(
      w,
      minN = minN,
      binSize = binSize,
      distrib = "betabinomial",
      b = k
    )
    
    sse_bbin <- sum(w_grad * (empirical - e_combined_sorted_binned[,2])^2)
    b_and_sse[counter + 1, 1] <- k
    b_and_sse[counter + 1, 2] <- sse_bbin
    labels[counter] = paste(
      "betabin,b=",
      signif(k, 2),
      "; SSE=",
      signif(sse_bbin, 2)
    )
    
    if (sse_bbin < sse || k == r_sta) {
      b_choice <- k
      sse <- sse_bbin
      e_combined_sorted_binned_cached <- e_combined_sorted_binned
      b_and_sse_cached <- b_and_sse
    } else if (sse_bbin > sse) {
      break
    }
    
    counter <- counter + 1
  }
  list(
    e_combined_sorted_binned = e_combined_sorted_binned_cached,
    b_and_sse = b_and_sse_cached,
    b_choice = b_choice,
    sse = sse,
    labels = labels,
    counter = counter
  )
}

#' @title Optimize the overdispersion parameter
#'
#' @description Minimize sse for betabinomials
#'
#' @param w_grad w_grad
#' @param b_and_sse b_and_sse
#' @param b_choice b_choice
#' @param empirical empirical
#' @param counter counter
#' @param minN minimum coverage level
#' @param p binomial probability of success parameter
#' @param binSize approximate number of bins
#' @param r_by r_by
#' @return list
#' @export
optimize_overdispersion_parameter <- function(
  w_grad,
  b_and_sse,
  b_choice,
  sse,
  empirical,
  counter,
  minN = 6,
  p = 0.5,
  binSize = 40,
  r_by = 0.1
) {
  flag <- 3
  if (b_choice >= 0.9) {
    flag <- 0
    newctr <- counter
  }
  while (flag) {
    r_sta <- max(0, b_choice - r_by)
    r_end <- b_choice + r_by
    r_by <- r_by / 2
    b_range <- seq(r_sta, r_end, by = r_by)
    labels <- matrix(0, nrow = 50, ncol = 1)
    newctr <- 1
    sse <- b_and_sse[1, 2]
    b_choice <- 0

    for (k in b_range) {
      e_combined_sorted_binned <- nulldistrib(
        w,
        minN = minN,
        p = p,
        binSize = binSize,
        distrib = "betabinomial",
        b = k
      )

      ## minimize sse for betabinomials
      sse_bbin <- sum(w_grad * (empirical - e_combined_sorted_binned[,2])^2)
      b_and_sse[(counter + 2), 1] <- k
      b_and_sse[(counter + 2), 2] <- sse_bbin
      labels[newctr] = paste(
        "betabin,b=",
        signif(k, 3),
        "; SSE=",
        signif(sse_bbin, 3)
      )
      
      if (sse_bbin < sse) {
        sse <- sse_bbin
        b_choice <- k 
      } else if (sse_bbin > sse) {
        break
      }
      
      counter <- counter + 1
      newctr <- newctr + 1
    }
    flag <- flag - 1
    labels = labels[1:(newctr + 1),]
    if (signif(b_and_sse[counter + 2, 2], 3) == signif(b_and_sse[counter + 1, 2], 3)) {
      flag <- 0 
    }
  }
  list(
    e_combined_sorted_binned = e_combined_sorted_binned,
    counter = counter
  )
}
