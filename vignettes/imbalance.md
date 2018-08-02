---
title: "Chen 2016 allelic imbalance"
author: "Anthony Aylward"
date: "2018-08-02"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{Vignette Title}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---



First, import the libraries we will use


```r
library(VGAM)
#> Loading required package: methods
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
accb[["p_betabin"]] <- p_betabin
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
fp_betabinomial <- data.frame(
  pval = p_thresh,
  FP_betabin = sapply(p_thresh, function(x) fp(w, p, x, "betabinomial", b))
)
head(fp_betabinomial)
#>    pval FP_betabin
#> 1 0.000     0.0000
#> 2 0.001   156.5733
#> 3 0.002   310.6079
#> 4 0.003   492.6134
#> 5 0.004   670.8305
#> 6 0.005   886.3256
```

FDR.txt


```r
tp_bin <- sapply(p_thresh, cutoff, y = p_bin) + 1
tp_betabin <- sapply(p_thresh, cutoff , y = p_betabin) + 1
fdr_bin <- fp_binomial[,2] / tp_bin
fdr_betabin <- fp_betabinomial[,2] / tp_betabin
p_choice_bin <- max(p_thresh[fdr_bin <= FDR_thresh])
p_choice_betabin <- max(p_thresh[fdr_betabin <= FDR_thresh])

fdr_choice_bin <- max(fdr_bin[fdr_bin <= FDR_thresh])
fdr_choice_betabin <- max(fdr_betabin[fdr_betabin <= FDR_thresh])
```

Bisection method to find p-value


```r
p_choice_bin_1 <- as.data.frame(
  bisect(
    p_bin,
    fp_bin[,2],
    p_choice_bin,
    fdr_choice_bin,
    FDR_thresh,
    step,
    "binomial",
    b = 0,
    w,
    p
  )
)
p_choice_betabin_1 <- as.data.frame(
  bisect(
    p_betabin,
    fp_betabinomial[,2],
    p_choice_betabin,
    fdr_choice_betabin,
    FDR_thresh,
    step,
    "betabinomial",
    b = b,
    w,
    p
  )
)
head(p_choice_betabin_1)
#>         V1         V2            V3
#> 1 0.004000 0.04903015 10.0000000000
#> 2 0.003950 0.04882447  0.0011755293
#> 3 0.003975 0.04888301  0.0011169918
#> 4 0.004000 0.04903373  0.0009662665
#> 5 0.004025 0.04934589  0.0006541134
#> 6 0.004050 0.05088008 -0.0008800844
```


```r
p_choice_bin_1 <- p_choice_bin_1[p_choice_bin_1[,3] > 0,]
p_choice_bin_2 <- p_choice_bin_1[nrow(p_choice_bin_1), 1]
p_choice_betabin_1 <- p_choice_betabin_1[p_choice_betabin_1[,3] > 0,]
p_choice_betabin_2 <- p_choice_betabin_1[nrow(p_choice_betabin_1), 1]
```

Formatting FDR_txt


```r
FDR_txt <- data.frame(
  pval = p_thresh,
  P_bin = tp_bin,
  FP_bin = fp_binomial[,2],
  FDR_bin = fdr_bin,
  P_betabin = tp_betabin,
  FP_betabin = fp_betabinomial[,2],
  FDR_betabin = fdr_betabin
)
FDR_txt = FDR_txt[-nrow(FDR_txt),]
head(FDR_txt)
#>    pval P_bin   FP_bin     FDR_bin P_betabin FP_betabin FDR_betabin
#> 1 0.000     6   0.0000 0.000000000         1     0.0000  0.00000000
#> 2 0.001 19784 160.4480 0.008109987     11112   156.5733  0.01409047
#> 3 0.002 22108 316.5243 0.014317185     12203   310.6079  0.02545341
#> 4 0.003 23789 482.0072 0.020261768     13109   492.6134  0.03757826
#> 5 0.004 25108 630.4576 0.025109831     13682   670.8305  0.04903015
#> 6 0.005 26464 824.1110 0.031140834     14456   886.3256  0.06131195
```

Take in counts.txt and filter by p.betabin


```r
interestingHets_betabinom = accb[accb[["p_betabin"]] <= p_choice_betabin,]
head(interestingHets_betabinom)
#>     chr  start    end     TF_indiv_accB ref alt cA cC cG cT   p.binomial
#> 1  chr1  91604  91605  SA1_NA19099_accB   C   T  0 45  0  0 5.684342e-14
#> 6  chr1 714018 714019 POL2_NA18505_accB   A   G  2  0 19  0 2.212524e-04
#> 7  chr1 714018 714019 POL2_NA19099_accB   A   G  2  0 19  0 2.212524e-04
#> 8  chr1 714018 714019  SA1_NA18505_accB   A   G 54  0 16  0 5.853956e-06
#> 11 chr1 762484 762485 RPB2_NA11894_accB   C   A  0 11  0  0 9.765625e-04
#> 15 chr1 762600 762601 RPB2_NA11894_accB   T   C  0 12  0  0 4.882812e-04
#>    p.betabinomial    p_betabin
#> 1    0.0003760213 2.173231e-09
#> 6    0.0031147471 1.475486e-03
#> 7    0.0025591824 1.475486e-03
#> 8    0.2424655815 2.451691e-03
#> 11   0.0045665860 2.334527e-03
#> 15   0.0029682809 1.368218e-03
```
