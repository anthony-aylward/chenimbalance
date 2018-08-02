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
FDR.thresh <- 0.05
p <- 0.5
b <- 0.01855469
```

Prepare the data


```r
cACGT <- accb[c("cA", "cC", "cG", "cT")]
lower <- apply(cACGT, 1, function(x) sort(x, partial = 3)[3])
higher <- apply(cACGT, 1, max)
total <- higher + lower
p.bin = apply(
  data.frame(2 * mapply(pbinom, lower, total, p)),
  1,
  function(x) min(x, 1)
)
p.betabin = apply(
  data.frame(
    2 * mapply(pbetabinom, lower, total, p, b)),
    1,
    function(x) min(x, 1)
)
head(p.bin)
#> [1] 5.684342e-14 6.029232e-01 6.397694e-01 1.000000e+00 1.000000e+00
#> [6] 2.212524e-04
```

Simulations


```r
step = 0.0001
p_thresh = data.frame(
  c(
    seq(0, 0.01, by = 0.001),
    seq(0.01, 0.1, by = 0.01)[-1],
    seq(0.1, 1, by = 0.1)[-1]
  )
)
```

Table of empirical counts


```r
w = as.data.frame(table(total), stringsAsFactors = F)
w = w[as.numeric(w[,1]) >= 6,]
```

FP


```r
fp_binomial = data.frame(
  pval = p_thresh,
  FP_Bin = apply(p_thresh, 1, function(x) fp(w, p, x, "binomial"))
)
head(fp_binomial)
#>   c.seq.0..0.01..by...0.001...seq.0.01..0.1..by...0.01...1...seq.0.1..
#> 1                                                                0.000
#> 2                                                                0.001
#> 3                                                                0.002
#> 4                                                                0.003
#> 5                                                                0.004
#> 6                                                                0.005
#>     FP_Bin
#> 1   0.0000
#> 2 160.4480
#> 3 316.5243
#> 4 482.0072
#> 5 630.4576
#> 6 824.1110
```
