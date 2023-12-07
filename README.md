# gedepir

**G**ene **E**xpression **DE**convolution **PI**peline in **R**

Package {gedepir} provides practical implementation of recommended guidelines for inference of cell-type proportions from reference-free deconvolution algorithms.

## Installation

To build a [conda](%60) environment :

``` bash
conda create -n rtmp -c conda-forge r-base r-biocmanager r-devtools r-ggplot2 r-rmarkdown r-dplyr r-cluster r-data.table r-fastica r-clue r-pheatmap 
```

To install dependencies and get the current version from GitHub:

``` r
BiocManager::install(c("NMF","fgsea"))
remotes::install_github("bcm-uga/gedepir")
```

## Usage

Vignette for basic usage is available at : <https://bcm-uga.github.io/gedepir/>
