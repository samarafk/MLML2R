---
title: "`MLML2R`: an R package for maximum likelihood estimates of DNA methylation and hydroxymethylation"
author: Samara F. Kiihl, Maria Tellez-Plaza
output:
  html_document:
    keep_md: yes
    toc: true
    number_sections: true
---



# Introduction


This document presents an example of the usage of the `MLML2R` package for R.

Install the R package using the following commands on the R console:


```r
install.packages("devtools")
devtools::install_github("samarafk/MLML2R")
library(MLML2R)
```


The function `MLML` provides maximum likelihood estimates (MLE) for 5-hmC and 5-mC levels using data from any combination of two of the methods: BS-seq, TAB-seq or oxBS-seq, or combination of all the three methods.

The algorithm implemented in the `MLML` function is based on the Expectation-Maximization (EM) algorithm proposed by [Qu *et al.* (2013)](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3789553/). In addition, our implementation is optimized, since we derived the exact constrained MLE for 5-mC or 5-hmC levels, and the iterative EM algorithm is not needed. Our improved formulation can, thus, decrease analytic processing time and computational burden, common bottlenecks when processing single-base profiling data from thousands of samples.

Furthermore, our routine is flexible and can be used with both next generation sequencing and Infinium Methylation microarray data in the R-statistical language.



# True proportion



```r
library(MLML2R)
```


True proportions used in the data simulation (`true_parameters_sim2` dataset).

<img src="README_files/figure-html/unnamed-chunk-3-1.png" style="display: block; margin: auto;" />


Dataset used in this example was simulated from the above true proportions.


# BS+oxBS+TAB data




```r
# via EM algortihm
results_em <- MLML(T = MethylatedBS_sim2 , U = UnMethylatedBS_sim2,
                   L = UnMethylatedOxBS_sim2, M = MethylatedOxBS_sim2,
                   G = UnMethylatedTAB_sim2, H = MethylatedTAB_sim2,tol=0.00001,iterative = TRUE)

# via Lagrange multiplier approximation
results_lag <- MLML(T = MethylatedBS_sim2 , U = UnMethylatedBS_sim2,
                   L = UnMethylatedOxBS_sim2, M = MethylatedOxBS_sim2,
                   G = UnMethylatedTAB_sim2, H = MethylatedTAB_sim2)
```

<img src="README_files/figure-html/unnamed-chunk-5-1.png" style="display: block; margin: auto;" />


```r
all.equal(results_em$hmC,results_lag$hmC)
```

```
## [1] "Mean relative difference: 1.371669e-05"
```

```r
all.equal(results_em$C,results_lag$C)
```

```
## [1] "Mean relative difference: 2.519594e-06"
```

```r
all.equal(results_em$hmC[,1],true_parameters_sim2$p_h)
```

```
## [1] "names for target but not for current"
## [2] "Mean relative difference: 0.04164895"
```

```r
all.equal(results_em$mC[,1],true_parameters_sim2$p_m)
```

```
## [1] "names for target but not for current"
## [2] "Mean relative difference: 0.01735317"
```

```r
all.equal(results_lag$hmC[,1],true_parameters_sim2$p_h)
```

```
## [1] "names for target but not for current"
## [2] "Mean relative difference: 0.04165114"
```

```r
all.equal(results_lag$mC[,1],true_parameters_sim2$p_m)
```

```
## [1] "names for target but not for current"
## [2] "Mean relative difference: 0.01737099"
```


```r
library(microbenchmark)
mbm = microbenchmark(
   lagrange = MLML(T = MethylatedBS_sim2 , U = UnMethylatedBS_sim2,
                                                             L = UnMethylatedOxBS_sim2, M = MethylatedOxBS_sim2,
                                                            G = UnMethylatedTAB_sim2, H = MethylatedTAB_sim2),
   EM = MLML(T = MethylatedBS_sim2 , U = UnMethylatedBS_sim2,
                            L = UnMethylatedOxBS_sim2, M = MethylatedOxBS_sim2,
                            G = UnMethylatedTAB_sim2, H = MethylatedTAB_sim2,tol=0.0001,iterative = TRUE),
   times=10)
mbm
```

```
## Unit: milliseconds
##      expr        min         lq      mean    median        uq       max
##  lagrange   8.797352   9.166617  21.61562  10.12498  49.17773  49.77772
##        EM 232.395691 240.551371 275.79024 276.31503 313.62644 313.86383
##  neval
##     10
##     10
```

```r
T = MethylatedBS_sim2
U = UnMethylatedBS_sim2
L = UnMethylatedOxBS_sim2
M = MethylatedOxBS_sim2
G = UnMethylatedTAB_sim2
H = MethylatedTAB_sim2

results_oxBS_TAB_BS_exact <- list()

results_oxBS_TAB_BS_exact$mC <- M/(M+L)
results_oxBS_TAB_BS_exact$hmC <- H/(H+G)
results_oxBS_TAB_BS_exact$C <- U/(U+T)


range(results_oxBS_TAB_BS_exact$mC)
```

```
## [1] 0.009537784 0.969560706
```

```r
range(results_oxBS_TAB_BS_exact$hmC)
```

```
## [1] 0.0000000 0.7659855
```

```r
range(results_oxBS_TAB_BS_exact$C)
```

```
## [1] 0.02216867 0.98640249
```

<img src="README_files/figure-html/unnamed-chunk-8-1.png" style="display: block; margin: auto;" />



<img src="README_files/figure-html/unnamed-chunk-9-1.png" style="display: block; margin: auto;" />


<img src="README_files/figure-html/unnamed-chunk-10-1.png" style="display: block; margin: auto;" />




# BS+oxBS data


```r
# obtain MLE via EM-algorithm for BS+oxBS:
results_em <- MLML(T = MethylatedBS_sim2 , U = UnMethylatedBS_sim2,
L = UnMethylatedOxBS_sim2, M = MethylatedOxBS_sim2,iterative=TRUE,tol=0.0001)

# obtain constrained exact MLE for BS+oxBS:
results_exact <- MLML(T = MethylatedBS_sim2 , U = UnMethylatedBS_sim2,
L = UnMethylatedOxBS_sim2, M = MethylatedOxBS_sim2)

all.equal(results_em$hmC,results_exact$hmC)
```

```
## [1] "Mean relative difference: 0.001189797"
```

<img src="README_files/figure-html/unnamed-chunk-12-1.png" style="display: block; margin: auto;" />


<img src="README_files/figure-html/unnamed-chunk-13-1.png" style="display: block; margin: auto;" />


# BS+TAB data



```r
# obtain MLE via EM-algorithm for BS+TAB:
results_em <- MLML(T = MethylatedBS_sim2 , U = UnMethylatedBS_sim2,
G = UnMethylatedTAB_sim2, H = MethylatedTAB_sim2,tol=0.0001)

# obtain constrained exact MLE for BS+TAB:
results_exact <- MLML(T = MethylatedBS_sim2 , U = UnMethylatedBS_sim2,
G = UnMethylatedTAB_sim2, H = MethylatedTAB_sim2)

all.equal(results_em$hmC,results_exact$hmC)
```

```
## [1] TRUE
```



<img src="README_files/figure-html/unnamed-chunk-15-1.png" style="display: block; margin: auto;" />


<img src="README_files/figure-html/unnamed-chunk-16-1.png" style="display: block; margin: auto;" />



# TAB+oxBS data



```r
# obtain MLE via EM-algorithm for oxBS+TAB:
results_em <- MLML(L = UnMethylatedOxBS_sim2, M = MethylatedOxBS_sim2,
G = UnMethylatedTAB_sim2, H = MethylatedTAB_sim2,iterative=TRUE,tol=0.0001)

# obtain constrained exact MLE for oxBS+TAB:
results_exact <- MLML(L = UnMethylatedOxBS_sim2, M = MethylatedOxBS_sim2,
G = UnMethylatedTAB_sim2, H = MethylatedTAB_sim2)

all.equal(results_em$hmC,results_exact$hmC)
```

```
## [1] "Attributes: < Component \"dimnames\": Component 2: target is NULL, current is character >"
## [2] "Mean relative difference: 5.673549e-06"
```



<img src="README_files/figure-html/unnamed-chunk-18-1.png" style="display: block; margin: auto;" />


<img src="README_files/figure-html/unnamed-chunk-19-1.png" style="display: block; margin: auto;" />



