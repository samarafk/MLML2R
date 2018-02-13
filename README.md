# `MLML2R`: Maximum Likelihood Estimation of DNA Methylation and Hydroxymethylation Levels

The function `MLML` provides maximum likelihood estimates (MLE) for 5-hmC and 5-mC levels using data from the methods: BS-seq, TAB-seq or oxBS-seq (all three methods or any combination of two methods). 

Estimates can be obtained using the Expectation-Maximization (EM) algorithm proposed by [Qu *et al.* (2013)](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3789553/), the pool-adjacent-violators algorithm (PAVA) or an approximated solution of the Lagrange multiplier method. 

Furthermore, our routine is flexible and can be used with both next generation sequencing and Infinium Methylation microarray data in the R-statistical language.


The package is available on CRAN:

[https://CRAN.R-project.org/package=MLML2R](https://CRAN.R-project.org/package=MLML2R).


To install and load the package, use the following commands on the R console:

```r
install.packages("MLML2R")
library(MLML2R)
```
