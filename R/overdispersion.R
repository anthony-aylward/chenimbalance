#===============================================================================
# overdispersion.R
#===============================================================================

# Imports ======================================================================

#' @import parallel




# Functions ====================================================================

#' @title Choose overdispersion parameter
#'
#' @description Minimize sse for betabinomials
#'
#' @param w_grad w_grad
#' @param w w
#' @param empirical empirical
#' @param sse sse
#' @param minN minimum coverage level
#' @param binSize approximate number of bins
#' @param r_sta,r_end start and end of the parameter range
#' @param r_by step of the parameter range
#' @param cores number of cores to use
#' @return list
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
  r_by = 0.1,
  p = 0.5,
  cores = detectCores()
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
  
  break_signal <- FALSE
  for (i in seq(to = length(b_range), by = cores)) {
    distribution_list <- mclapply(
      b_range[i:min(length(b_range), i + cores - 1)],
      function(k) {
        nulldistrib(
          w,
          minN = minN,
          binSize = binSize,
          distrib = "betabinomial",
          p = p,
          b = k
        )
      },
      mc.cores = cores
    )
    for (j in 1:min(cores, length(b_range) - i)) {
      k <- b_range[[i + j - 1]]
      e_combined_sorted_binned <- distribution_list[[j]]
      
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
      } else if (sse_bbin > sse) {
        break_signal = TRUE
        break
      }

      counter <- counter + 1
    }
    if (break_signal) {
      break
    }
  }
  list(
    e_combined_sorted_binned = e_combined_sorted_binned_cached,
    b_and_sse = b_and_sse,
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
#' @param w w
#' @param b_and_sse b_and_sse
#' @param b_choice b_choice
#' @param empirical empirical
#' @param counter counter
#' @param minN minimum coverage level
#' @param p binomial probability of success parameter
#' @param binSize approximate number of bins
#' @param r_by r_by
#' @param cores number of cores to use
#' @return list
#' @export
optimize_overdispersion_parameter <- function(
  w_grad,
  w,
  b_and_sse,
  b_choice,
  sse,
  empirical,
  counter,
  minN = 6,
  p = 0.5,
  binSize = 40,
  r_by = 0.1,
  cores = detectCores()
) {
  flag <- 3
  if (b_choice >= 0.9) {
    flag <- FALSE
    newctr <- counter
  }
  while (flag > 0) {
    r_sta <- max(0, b_choice - r_by)
    r_end <- b_choice + r_by
    r_by <- r_by / 2
    b_range <- seq(r_sta, r_end, by = r_by)
    labels <- matrix(0, nrow = 50, ncol = 1)
    newctr <- 1
    
    break_signal <- FALSE
    for (i in seq(to = length(b_range), by = cores)) {
      distribution_list <- mclapply(
        b_range[i:min(length(b_range), i + cores - 1)],
        function(k) {
          nulldistrib(
            w,
            minN = minN,
            binSize = binSize,
            distrib = "betabinomial",
            p = p,
            b = k
          )
        },
        mc.cores = cores
      )
      for (j in 1:min(cores, length(b_range) - i)) {
        k <- b_range[[i + j - 1]]
        e_combined_sorted_binned <- distribution_list[[j]]
        sse_bbin <- sum(w_grad * (empirical - e_combined_sorted_binned[,2])^2)

        if (sse_bbin < sse) {
          sse <- sse_bbin
          b_choice <- k 
        }
      }
    }
    b_and_sse[(counter + 2), 1] <- b_choice
    b_and_sse[(counter + 2), 2] <- sse
    labels[newctr] = paste(
      "betabin,b=",
      signif(k, 3),
      "; SSE=",
      signif(sse_bbin, 3)
    )
    labels = labels[1:(newctr + 1),]
    if (
      signif(b_and_sse[counter + 2, 2], 3)
      == signif(b_and_sse[counter + 1, 2], 3)
    ) {
      flag <- flag - 1
    }
    counter <- counter + 1
    newctr <- newctr + 1
  }
  list(
    e_combined_sorted_binned = e_combined_sorted_binned,
    b_and_sse = b_and_sse,
    b_choice = b_choice,
    sse = sse,
    counter = counter
  )
}
