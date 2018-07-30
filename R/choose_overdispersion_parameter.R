#===============================================================================
# choose_overdispersion_parameter.R
#===============================================================================

#' @title Empirical allelic Ratio
#'
#' @details
#' Plot right closed interval left open (range] for x axis
#' note that you have pseudozeroes as counts in your data so thats fine
#'
#' @param data A data frame with columns "total" and "allelicRatio"
#' @param bins Breakpoints for bins
#' @return Histogram of allelic ratios
#' @export
# choose_overdispersion_parameter <- function(
#   w,
#   minN = 6,
#   binSize = 40,
# ) {
#   w <- weight_by_empirical_counts(total)
#   for (k in b.range) {
#     e.combined.sorted.binned = nulldistrib(
#       w,
#       minN = minN
#       binSize,
#       distrib = "betabinomial",
#       b = k
#     )
#   #   par(new=TRUE)
#   #   plot(e.combined.sorted.binned,ylim=c(0,yuplimit),pch=16,type='b',col=colors[ctr],bty='n',ylab='',xlab='',yaxt='n',xaxt='n',yaxs="i")
    
#     ## minimize sse for betabinomials
#     if (b.choice == 0) {
#       b.and.sse[1, 1] <- b.choice
#       b.and.sse[1, 2] <- sse
#     }
#     sse.bbin <- sum(w.grad * ((empirical-e.combined.sorted.binned[,2])^2))
#     b.and.sse[ctr+1,1] = k
#     b.and.sse[ctr+1,2] = sse.bbin
#     labels[ctr] = paste("betabin,b=",signif(k,2),"; SSE=",signif(sse.bbin,2))
    
#     if(sse.bbin < sse){ sse = sse.bbin; b.choice = k }
#     else if(sse.bbin > sse){ break }
    
#     ctr = ctr + 1
#   }
# }