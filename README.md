# scMNMF (single-cell modular Non-negative Matrix Factorization)

![R](https://img.shields.io/badge/R-%3E%3D4.0.0-blue.svg)
![License](https://img.shields.io/badge/license-GPL--3-green.svg)

`scMNMF` is an R package designed for single-cell RNA-seq data integration and interpretable biological signal decoupling. By incorporating graph-regularized constraints to explicitly model multiple biological variables (such as batches, biological conditions, and cell types), scMNMF effectively disentangles common biological programs from condition-specific variations.

## Installation

To run `scMNMF`, you can install it directly from GitHub using the `devtools` package:

```R
# Install devtools if not already installed
if (!requireNamespace("devtools", quietly = TRUE))
    install.packages("devtools")

# Install scMNMF from GitHub
library(devtools)
devtools::install_github("zhxran/scMNMF")
```

## **Dependencies**

**scMNMF** has been successfully tested on Windows, Linux, and Mac OS (R version >= 4.0.0). The package relies on the following dependencies:

* **Core Logic:** `Rcpp`, `RcppEigen` (for high-performance matrix operations)
* **Numerical Computing:** `Matrix`, `irlba`, `RSpectra`, `stats`, `base`
* **Data Manipulation:** `dplyr`, `matrixStats`, `sparseMatrixStats`
* **Bioinformatics Tools:** `Seurat` (>= 4.0.0), `mclust`, `BiocNeighbors`
* **Visualization:** `ggplot2`, `patchwork`, `pheatmap`
