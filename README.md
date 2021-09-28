# gedepir
Gene Expression DEconvolution PIpeline in R

Package {gedepir} (Gene Expression DEconvoluation PIpeline in R) provides practical implementation of recommended guidelines for inference of cell-type proportions from reference-free deconvolution algorithms. 

## Installation
To build a [conda](`) environment :
```bash
conda create -n rtmp -c conda-forge r-base r-biocmanager r-devtools r-ggplot2 r-rmarkdown r-dplyr r-cluster r-data.table r-fastica r-clue r-pheatmap 

```
To install dependencies and get the current version from GitHub:

```R
BiocManager::install(c("NMF","fgsea"))
remotes::install_github("bcm-uga/gedepir",build_vignettes = TRUE)
```


