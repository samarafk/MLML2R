# `MLML2R`: An R Package for Maximum Likelihood Estimation of DNA Methylation and Hydroxymethylation Proportions

Samara F. Kiihl, Maria Tellez-Plaza  


This document presents an example of the usage of the `MLML2R` package for R.


Proposed analyses of single-base profiling of either 5-hmC or 5-mC require combining data obtained using bisulfite conversion, oxidative bisulfite conversion or Tet-Assisted bisulfite conversion methods, but doing so naively produces inconsistent estimates of 5-mC or 5-hmC level [(Qu *et al.*, 2013)](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3789553/). 


The function `MLML` provides maximum likelihood estimates (MLE) for 5-hmC and 5-mC levels using data from any combination of two of the methods: BS-seq, TAB-seq or oxBS-seq, or combination of all the three methods.

The algorithm implemented in the `MLML` function is based on the Expectation-Maximization (EM) algorithm proposed by [Qu *et al.* (2013)](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3789553/). In addition, our implementation is optimized, since we derived the exact constrained MLE for 5-mC or 5-hmC levels, and the iterative EM algorithm is not needed. Our improved formulation can, thus, decrease analytic processing time and computational burden, common bottlenecks when processing single-base profiling data from thousands of samples.

Furthermore, our routine is flexible and can be used with both next generation sequencing and Infinium Methylation microarray data in the R-statistical language.



To install the R package, use the following commands on the R console:


```r
install.packages("devtools")
devtools::install_github("samarafk/MLML2R")
library(MLML2R)
```

Examples of usage:

* [Example 1](./data-raw/example1/): dataset from [Field *et al.* (2015)](https://doi.org/10.1371/journal.pone.0118202), which consists of eight DNA samples from the same DNA source treated with oxBS-BS and hybridized to the Infinium 450K array.

* [Example 2](./data-raw/example2/): dataset from [Johnson *et al.* (2016)](https://www.ncbi.nlm.nih.gov/pubmed/27886174), which consists of 30 DNA samples from glioblastoma tumors treated with oxBS-BS and hybridized to the Infinium 450K array.

* [Example 3](./data-raw/example3/): simulated dataset included in the package, which consists of 4 samples from BS, oxBS and TAB.
