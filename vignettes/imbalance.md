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




```r
FDR.thresh <- 0.05
p <- 0.5
cACGT <- data.frame(cA = accb$cA, cC = accb$cC, cG = accb$cG, cT = accb$cT)
#> Error in data.frame(cA = accb$cA, cC = accb$cC, cG = accb$cG, cT = accb$cT): object 'accb' not found
lower <- apply(cACGT, 1, function(x) sort(x, partial = 3)[3])
#> Error in apply(cACGT, 1, function(x) sort(x, partial = 3)[3]): object 'cACGT' not found
higher <- apply(cACGT, 1, max)
#> Error in apply(cACGT, 1, max): object 'cACGT' not found
total <- higher + lower
#> Error in eval(expr, envir, enclos): object 'higher' not found
b <- 0.01855469
p.bin = apply(
  data.frame(2 * mapply(pbinom, lower, total, p)),
  1,
  function(x) min(x, 1)
)
#> Error in mapply(pbinom, lower, total, p): object 'lower' not found
head(p.bin)
#> Error in head(p.bin): object 'p.bin' not found
```
