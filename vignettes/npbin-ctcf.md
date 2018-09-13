---
title: "Estimate overdispersion on the ctcf dataset from the NPBin paper"
author: "Anthony Aylward"
date: "2018-09-13"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{Vignette Title}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---



## Prepare the data.

This vignette requires the [npbin](https://github.com/anthony-aylward/npbin) 
package.


```r
library(chenimbalance)
library(npbin)
total_reads <- ctcf[["m"]]
data <- data.frame(
  total = total_reads,
  allelicRatio =  ctcf[["xm"]] / ctcf[["m"]]
)
data <- data[total_reads >= 5,]
data <- data[1:2000,]
head(data)
#>   total allelicRatio
#> 1   190    0.5368421
#> 2   185    0.4864865
#> 3    16    0.5625000
#> 4   375    0.2640000
#> 5    22    1.0000000
#> 6    40    0.6750000
```

## Empirical distribution

Compute the empirical allelic ratio distribution


```r
binSize <- 40
bins <- pretty(0:1, binSize)
minN <- 6
maxN <- min(2500, max(data[["total"]]))
empirical <- empirical_allelic_ratio(
  data,
  bins,
  maxN = maxN,
  minN = minN,
  plot = TRUE
)
```

![plot of chunk ctcf_empirical_dist](figure/ctcf_empirical_dist-1.png)

## Expected Binomial and Beta-Binomial distributions

Compute the weighted expected binomial distribution


```r
w <- weight_by_empirical_counts(data[["total"]])
d_combined_sorted_binned <- nulldistrib(
  w,
  minN = minN,
  binSize = binSize
)
```

Compute the sum of squared errors for the empirical distribution vs the 
weighted expected binomial distribution.


```r
sse <- sum((empirical - d_combined_sorted_binned[,2])^2)
sse
#> [1] 0.007432315
```

Choose the overdispersion parameter for the beta-binomial distribution


```r
w_grad <- graded_weights_for_sse_calculation(r_min = 0, r_max = 1, bins = bins)
overdispersion_details <- choose_overdispersion_parameter(
  w_grad,
  w,
  empirical,
  sse
)
head(overdispersion_details[["b_and_sse"]])
#>        b         sse
#> [1,] 0.0 0.007432315
#> [2,] 0.1 0.015246404
#> [3,] 0.2 0.024126144
#> [4,] 0.0 0.000000000
#> [5,] 0.0 0.000000000
#> [6,] 0.0 0.000000000
```

Generate a plot of the weighted expected binomial and weighted expected 
beta-binomial distributions overlaid on the empirical distribution


```r
plot_distributions(
  minN,
  maxN,
  bins,
  empirical,
  d_combined_sorted_binned,
  overdispersion_details[["e_combined_sorted_binned"]],
  yuplimit = 0.15
)
```

![plot of chunk ctcf_plot_initial_dist](figure/ctcf_plot_initial_dist-1.png)

`overdispersion_details` is a list whose elements include the chosen value of 
`b` and the sum of squared errors.


```r
paste(
  "b_chosen =",
  overdispersion_details[["b_choice"]],
  ", SSE_chosen =",
  overdispersion_details[["sse"]]
)
#> [1] "b_chosen = 0.1 , SSE_chosen = 0.0152464044437993"
```

Optimize the overdispersion parameter


```r
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
plot_distributions(
  minN,
  maxN,
  bins,
  empirical,
  d_combined_sorted_binned,
  optimized_overdispersion_details[["e_combined_sorted_binned"]],
  yuplimit = 0.15
)
```

![plot of chunk ctcf_optimize_parameter](figure/ctcf_optimize_parameter-1.png)

Check the optimized value


```r
list(
  b = optimized_overdispersion_details[["b_choice"]],
  sse = optimized_overdispersion_details[["sse"]]
)
#> $b
#> [1] 0.008984375
#> 
#> $sse
#> [1] 0.001672135
```

Plot the parameter search space


```r
b_and_sse <- (
  optimized_overdispersion_details[["b_and_sse"]]
  [1:(optimized_overdispersion_details[["counter"]] + 2),]
)
plot(
  b_and_sse[order(b_and_sse[,1]),],
  type = "b",
  pch = 16,
  xlim = c(min(b_and_sse[,1]), max(b_and_sse[,1])),
  ylim = c(min(b_and_sse[,2]), max(b_and_sse[,2]))
)
par(new = TRUE)
plot(
  optimized_overdispersion_details[["b_choice"]],
  optimized_overdispersion_details[["sse"]],
  bty = "n",
  ylab = "",
  xlab = "",
  yaxt = "n",
  xaxt = "n",
  col = "red",
  pch = 8,
  xlim = c(min(b_and_sse[,1]), max(b_and_sse[,1])),
  ylim = c(min(b_and_sse[,2]), max(b_and_sse[,2]))
)
```

![plot of chunk ctcf_plot_space](figure/ctcf_plot_space-1.png)

Compute the symmetric shape parameter and plot the estimated null beta (gold)
superimposed with the null beta estimated from NPBin (blue).


```r
shape = 1 / (2 * optimized_overdispersion_details[["b_choice"]]) - 1 / 2
plot_estimated_null(
  data[["allelicRatio"]],
  shape1_shape2 = c(37.53662, 37.30897),
  shape3_shape4 = c(shape, shape)
)
```

![plot of chunk npbin_ctcf_estimated_null](figure/npbin_ctcf_estimated_null-1.png)
