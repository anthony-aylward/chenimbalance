#===============================================================================
# choose_overdispersion_parameter.R
#===============================================================================

#' @title Provide a color palette
#'
#' @return A character vector representing the color palette
#' @export
#' @seealso \code{\link{choose_overdispersion_parameter}}
color_palette <- function() {
  c(
    "green", "blue", "orange","cyan", "pink", "purple", "brown", "black",
    "slategray1", "violetred", "tan", "deeppink", "darkgreen", "orchid",
    "darksalmon", "antiquewhite3", "magenta", "darkblue", "peru", "slateblue",
    "thistle", "tomato", "rosybrown1", "royalblue", "olivedrab"
  )
}

#' Choose overdispersion parameter
#'
#' Minimize sse for betabinomials
#'
#' @param data A data frame with columns "total" and "allelicRatio"
#' @param bins Breakpoints for bins
#' @return Histogram of allelic ratios
#' @export
choose_overdispersion_parameter <- function(
  w.grad,
  w,
  empirical,
  sse,
  minN = 6,
  binSize = 40,
  b_choice = 0,
  r_sta = 0,
  r_end = 0.99,
  r_by  = 0.1
) {
  colors <- color_palette()
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
  #   par(new=TRUE)
  #   plot(e.combined.sorted.binned,ylim=c(0,yuplimit),pch=16,type='b',col=colors[counter],bty='n',ylab='',xlab='',yaxt='n',xaxt='n',yaxs="i")
    
    sse_bbin <- sum(w_grad * (empirical - e_combined_sorted_binned[,2])^2)
    b_and_sse[counter + 1, 1] <- k
    b_and_sse[counter + 1, 2] <- sse_bbin
    labels[counter] = paste(
      "betabin,b=",
      signif(k, 2),
      "; SSE=",
      signif(sse_bbin, 2)
    )
    
    if (sse_bbin < sse) {
      b_choice <- k
      sse <- sse_bbin
    } else if (sse_bbin > sse) {
      break
    }
    
    counter <- counter + 1
  }
  list(
    e_combined_sorted_binned = e_combined_sorted_binned,
    b_and_sse = b_and_sse,
    b_choice = b_choice,
    sse = sse,
    labsls = labels
  )
}