---
title: "Chen 2016 allelic imbalance"
author: "Anthony Aylward"
date: "2018-08-01"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{Vignette Title}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---



First, import the libraries we will use


```r
library(VGAM)
#> Loading required package: stats4
#> Loading required package: splines
library(chenimbalance)
```

Set some parameters


```r
FDR_thresh <- 0.05
p <- 0.5
b <- 0.01855469
```

Prepare the data


```r
cACGT <- accb[c("cA", "cC", "cG", "cT")]
lower <- apply(cACGT, 1, function(x) sort(x, partial = 3)[3])
higher <- apply(cACGT, 1, max)
total <- higher + lower
p_bin <- apply(
  data.frame(2 * mapply(pbinom, lower, total, p)),
  1,
  function(x) min(x, 1)
)
p_betabin <- apply(
  data.frame(
    2 * mapply(pbetabinom, lower, total, p, b)),
    1,
    function(x) min(x, 1)
)
head(p_bin)
#> [1] 5.684342e-14 6.029232e-01 6.397694e-01 1.000000e+00 1.000000e+00
#> [6] 2.212524e-04
```

Simulations


```r
step <- 0.0001
p_thresh <- c(
  seq(0, 0.01, by = 0.001),
  seq(0.01, 0.1, by = 0.01)[-1],
  seq(0.1, 1, by = 0.1)[-1]
)
```

Table of empirical counts


```r
w <- as.data.frame(table(total), stringsAsFactors = F)
w <- w[as.numeric(w[,1]) >= 6,]
```

FP


```r
fp_binomial <- data.frame(
  pval = p_thresh,
  FP_bin = sapply(p_thresh, function(x) fp(w, p, x, "binomial"))
)
head(fp_binomial)
#>    pval   FP_bin
#> 1 0.000   0.0000
#> 2 0.001 160.4480
#> 3 0.002 316.5243
#> 4 0.003 482.0072
#> 5 0.004 630.4576
#> 6 0.005 824.1110
```


```r
fp_betabinomial <- as.data.frame(
  pval = p_thresh,
  FP_betabin = sapply(p_thresh, function(x) fp(w, p, x, "betabinomial", b))
)
#> Error in as.data.frame(pval = p_thresh, FP_betabin = sapply(p_thresh, : argument "x" is missing, with no default
head(fp_betabinomial)
#> Error in head(fp_betabinomial): object 'fp_betabinomial' not found
```
