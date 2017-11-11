---
title: "Example data for `MLML2R` package"
output:
  html_document:
    keep_md: yes
    toc: true
    number_sections: true
---





# Getting publicly available data


We will use the dataset from [Field *et al.* (2015)](https://doi.org/10.1371/journal.pone.0118202), which consists of eight DNA samples from the same DNA source treated with oxBS-BS and hybridized to the Infinium 450K array.

The steps shown in this section follows the [vignette](https://kasperdanielhansen.github.io/genbioconductor/html/minfi.html) from `minfi` package.

We start with the steps to get the raw data from the GEO repository.
The dataset from [Field *et al.* (2015)](https://doi.org/10.1371/journal.pone.0118202) is available at GEO accession [GSE63179](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE63179).

The sample was divided into four BS and four oxBS replicates.
 	
Platform used: GPL16304	Illumina HumanMethylation450 BeadChip [UBC enhanced annotation v1.0]

Samples:

* GSM1543269	brain-BS-1
* GSM1543270	brain-oxBS-3
* GSM1543271	brain-BS-2
* GSM1543272	brain-oxBS-4
* GSM1543273	brain-BS-3
* GSM1543274	brain-BS-4
* GSM1543275	brain-oxBS-1
* GSM1543276	brain-oxBS-2



This example has the following dependencies:


```r
library(minfi)
library(GEOquery)
```

Use the following commands to install these packages in R:

```r
source("http://www.bioconductor.org/biocLite.R")
biocLite(c("minfi", "GEOquery"))
```






```r
getGEOSuppFiles("GSE63179")
untar("GSE63179/GSE63179_RAW.tar", exdir = "GSE63179/idat")
head(list.files("GSE63179/idat", pattern = "idat"))
```


Decompress the compressed IDAT files:


```r
idatFiles <- list.files("GSE63179/idat", pattern = "idat.gz$", full = TRUE)
sapply(idatFiles, gunzip, overwrite = TRUE)
```

```
## named list()
```


Now we read the IDAT files in the directory:


```r
rgSet <- read.metharray.exp("GSE63179/idat")
rgSet
```

```
## class: RGChannelSet 
## dim: 622399 8 
## metadata(0):
## assays(2): Green Red
## rownames(622399): 10600313 10600322 ... 74810490 74810492
## rowData names(0):
## colnames(8): GSM1543269_9373551079_R01C01
##   GSM1543270_9373551079_R01C02 ... GSM1543275_9373551079_R05C01
##   GSM1543276_9373551079_R06C01
## colData names(0):
## Annotation
##   array: IlluminaHumanMethylation450k
##   annotation: ilmn12.hg19
```


```r
pData(rgSet)
```

```
## DataFrame with 8 rows and 0 columns
```


```r
sampleNames(rgSet)
```

```
## [1] "GSM1543269_9373551079_R01C01" "GSM1543270_9373551079_R01C02"
## [3] "GSM1543271_9373551079_R02C01" "GSM1543272_9373551079_R02C02"
## [5] "GSM1543273_9373551079_R03C01" "GSM1543274_9373551079_R04C01"
## [7] "GSM1543275_9373551079_R05C01" "GSM1543276_9373551079_R06C01"
```

The file names consists of a GEO identifier (the GSM part) followed by a standard IDAT naming convention with a 10 digit number which is an array identifier followed by an identifier of the form R01C01. This is because each array actually allows for the hybridization of 12 samples in a 6x2 arrangement. The 9373551079_R01C01 means row 1 and column 1 on chip 9373551079. 

We need to identify the samples from different methods: BS-conversion, oxBS-conversion.





```r
geoMat <- getGEO("GSE63179")
pD.all <- pData(geoMat[[1]])
pD <- pD.all[, c("title", "geo_accession", "characteristics_ch1.1", "characteristics_ch1.2","characteristics_ch1.3")]
pD
```


```
##                   title geo_accession characteristics_ch1.1
## GSM1543269   brain-BS-1    GSM1543269          gender: Male
## GSM1543270 brain-oxBS-3    GSM1543270          gender: Male
## GSM1543271   brain-BS-2    GSM1543271          gender: Male
## GSM1543272 brain-oxBS-4    GSM1543272          gender: Male
## GSM1543273   brain-BS-3    GSM1543273          gender: Male
## GSM1543274   brain-BS-4    GSM1543274          gender: Male
## GSM1543275 brain-oxBS-1    GSM1543275          gender: Male
## GSM1543276 brain-oxBS-2    GSM1543276          gender: Male
##            characteristics_ch1.2 characteristics_ch1.3
## GSM1543269               age: 25    bisulfite_proc: BS
## GSM1543270               age: 25  bisulfite_proc: oxBS
## GSM1543271               age: 25    bisulfite_proc: BS
## GSM1543272               age: 25  bisulfite_proc: oxBS
## GSM1543273               age: 25    bisulfite_proc: BS
## GSM1543274               age: 25    bisulfite_proc: BS
## GSM1543275               age: 25  bisulfite_proc: oxBS
## GSM1543276               age: 25  bisulfite_proc: oxBS
```


```r
names(pD)[c(3,4,5)] <- c("gender", "age","method")
pD$gender <- sub("^gender: ", "", pD$gender)
pD$age <- sub("^age: ", "", pD$age)
pD$method <- sub("^bisulfite_proc: ","",pD$method)
```

We now need to merge this pheno data into the methylation data. The following are commands to make sure we have the same row identifier in both datasets before merging.


```r
sampleNames(rgSet) <- sapply(sampleNames(rgSet),function(x) strsplit(x,"_")[[1]][1])
rownames(pD) <- pD$geo_accession
pD <- pD[sampleNames(rgSet),]
pData(rgSet) <- as(pD,"DataFrame")
rgSet
```

```
## class: RGChannelSet 
## dim: 622399 8 
## metadata(0):
## assays(2): Green Red
## rownames(622399): 10600313 10600322 ... 74810490 74810492
## rowData names(0):
## colnames(8): GSM1543269 GSM1543270 ... GSM1543275 GSM1543276
## colData names(5): title geo_accession gender age method
## Annotation
##   array: IlluminaHumanMethylation450k
##   annotation: ilmn12.hg19
```

```r
save(rgSet,file="rgSet.rds")
```


# Preprocessing

We refer the reader to the `minfi` package tutorials for more preprocessing options.


We need to install the required package bellow:

```r
source("https://bioconductor.org/biocLite.R")
biocLite("IlluminaHumanMethylation450kmanifest")
```


First, we removed probes with detection p-value <0.01 in any of the 8 arrays. The function `detectionP` identifies failed positions defined as both the methylated and unmethylated channel reporting background signal levels.


```r
detP <- detectionP(rgSet)
failed <- detP >0.01
## Keep probes which failed in at most maxFail arrays (0 = the probe passed in all arrays)
maxFail<- 0
keep_probes <- rowSums(failed) <= maxFail
```

We kept $83\%$ of the probes according to this criterion.


The `rgSet` object is a class called `RGChannelSet` which represents two color data with a green and a red channel. We will use, as input in the `MLML` funcion, a `MethylSet`, which contains the methylated and unmethylated signals. The most basic way to construct a `MethylSet` is to using the function `preprocessRaw` which uses the array design to match up the different probes and color channels to construct the methylated and unmethylated signals. Here we will use the `preprocessNoob` function, which does the preprocessing and returns a `MethylSet`.




Arrays were then normalized using the Noob/ssNoob preprocessing method for Infinium methylation microarrays.


From a `MethylSet` it is easy to compute Beta values, defined as:

Beta = Meth / (Meth + Unmeth + c)

The c constant is chosen to avoid dividing with small values. Illumina uses a default of c=100. The function `getBeta` from `minfi` package can be used to obtain the Beta values.



```r
MSet.noob<- preprocessNoob(rgSet[keep_probes,])
```

```
## [preprocessNoob] Applying R/G ratio flip to fix dye bias...
```





Prepare de input data:



```r
MethylatedBS <- getMeth(MSet.noob)[,c(1,3,5,6)]
UnMethylatedBS <- getUnmeth(MSet.noob)[,c(1,3,5,6)]

MethylatedOxBS <- getMeth(MSet.noob)[,c(7,8,2,4)]
UnMethylatedOxBS <- getUnmeth(MSet.noob)[,c(7,8,2,4)]
```

Only a small sample of CpGs and 2 replicates will be used in the example data.


```r
set.seed(112017)

a <- sample(1:dim(MethylatedBS)[1],100,replace=FALSE)

MethylatedBS <- MethylatedBS[a,1:2]
UnMethylatedBS <- UnMethylatedBS[a,1:2]
MethylatedOxBS <- MethylatedOxBS[a,1:2]
UnMethylatedOxBS <- UnMethylatedOxBS[a,1:2]
```




