#===============================================================================
# prob.R
#===============================================================================

# Imports ======================================================================

#' @import parallel




# Functions ====================================================================

#' @title Choose probability of success parameter
#'
#' @description Minimize sse for betabinomials
#'
#' @param w_grad w_grad
#' @param w w
#' @param empirical empirical
#' @param sse sse
#' @param b overdispersion parameter
#' @param minN minimum coverage level
#' @param binSize approximate number of bins
#' @param r_sta,r_end start and end of the parameter range
#' @param r_by step of the parameter range
#' @param n_cores number of cores to use
#' @return list
#' @export
choose_probability_of_success_parameter <- function(
  w_grad,
  w,
  empirical,
  sse,
  b,
  minN = 6,
  binSize = 40,
  prob_choice = 0.5,
  r_sta = 0.1,
  r_end = 0.99,
  r_by = 0.1,
  n_cores = detectCores()
) {
  counter <- 1
  prob_and_sse = matrix(
    c(prob_choice, sse, rep(0, 48)),
    nrow = 50,
    ncol = 2,
    byrow = TRUE,
    dimnames = list(NULL,  c("prob", "sse"))
  )
  labels <- matrix(0, nrow = 50, ncol = 1)
  prob_range <- seq(r_sta, r_end, by = r_by)
  
  break_signal <- FALSE
  for (i in seq(to = length(prob_range), by = n_cores)) {
    distribution_list <- mclapply(
      prob_range[i:min(length(prob_range), i + n_cores - 1)],
      function(k) {
        nulldistrib(
          w,
          minN = minN,
          binSize = binSize,
          distrib = "betabinomial",
          b = b,
          p = k
        )
      },
      mc.cores = n_cores
    )
    for (j in 1:min(n_cores, length(prob_range) - i)) {
      k <- prob_range[[i + j - 1]]
      e_combined_sorted_binned <- distribution_list[[j]]
      
      sse_bbin <- sum(w_grad * (empirical - e_combined_sorted_binned[,2])^2)
      prob_and_sse[counter + 1, 1] <- k
      prob_and_sse[counter + 1, 2] <- sse_bbin
      labels[counter] = paste(
        "betabin,prob=",
        signif(k, 2),
        "; SSE=",
        signif(sse_bbin, 2)
      )
      
      if (sse_bbin < sse || k == r_sta) {
        prob_choice <- k
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
    prob_and_sse = prob_and_sse,
    prob_choice = prob_choice,
    sse = sse,
    labels = labels,
    counter = counter
  )
}

#' @title Optimize the probability of success parameter
#'
#' @description Minimize sse for betabinomials
#'
#' @param w_grad w_grad
#' @param w w
#' @param prob_and_sse b_and_sse
#' @param prob_choice b_choice
#' @param empirical empirical
#' @param counter counter
#' @param b overdispersion parameter
#' @param minN minimum coverage level
#' @param binSize approximate number of bins
#' @param r_by r_by
#' @param n_cores number of cores to use
#' @return list
#' @export
optimize_probability_of_success_parameter <- function(
  w_grad,
  w,
  prob_and_sse,
  prob_choice,
  sse,
  empirical,
  counter,
  b,
  minN = 6,
  binSize = 40,
  r_by = 0.1,
  n_cores = detectCores()
) {
  flag <- TRUE
  if (prob_choice >= 0.9) {
    flag <- FALSE
    newctr <- counter
  }
  while (flag) {
    r_sta <- max(0, prob_choice - r_by)
    r_end <- prob_choice + r_by
    r_by <- r_by / 2
    prob_range <- seq(r_sta, r_end, by = r_by)
    labels <- matrix(0, nrow = 50, ncol = 1)
    newctr <- 1
    sse <- prob_and_sse[counter, 2]

    prob_and_sse[counter + 2,] <- matrix(c(prob_choice, sse), nrow = 1)
    counter <- counter + 1
    newctr <- newctr + 1
    
    break_signal <- FALSE
    for (i in seq(to = length(prob_range), by = n_cores)) {
      distribution_list <- mclapply(
        prob_range[i:min(length(prob_range), i + n_cores - 1)],
        function(k) {
          nulldistrib(
            w,
            minN = minN,
            binSize = binSize,
            distrib = "betabinomial",
            p = k,
            b = b
          )
        },
        mc.cores = n_cores
      )
      for (j in 1:min(n_cores, length(prob_range) - i)) {
        k <- prob_range[[i + j - 1]]
        e_combined_sorted_binned <- distribution_list[[j]]
        sse_bbin <- sum(w_grad * (empirical - e_combined_sorted_binned[,2])^2)

        if (sse_bbin < sse) {
          prob_and_sse[(counter + 2), 1] <- prob_choice
          prob_and_sse[(counter + 2), 2] <- sse
          labels[newctr] = paste(
            "betabin,prob=",
            signif(prob_choice, 3),
            "; SSE=",
            signif(sse, 3)
          )
          sse <- sse_bbin
          prob_choice <- k
          counter <- counter + 1
          newctr <- newctr + 1
        }
      }
    }
    labels = labels[1:(newctr + 1),]
    if (signif(prob_and_sse[counter + 2, 2], 3) == signif(prob_and_sse[counter + 1, 2], 3)) {
      flag <- FALSE
    }
  }
  list(
    e_combined_sorted_binned = e_combined_sorted_binned,
    prob_and_sse = prob_and_sse,
    prob_choice = prob_choice,
    sse = sse,
    counter = counter
  )
}
