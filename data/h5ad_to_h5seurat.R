#!/usr/bin/env Rscript

# Set CRAN mirror
options(repos = c(CRAN = "https://cloud.r-project.org"))

# Load required packages
install.packages("remotes")
remotes::install_github("mojaveazure/seurat-disk")
install.packages("hdf5r", dependencies = TRUE)

library(SeuratDisk, verbose = TRUE)
library(Seurat)


# Set working directory
setwd("/div/pythagoras/u1/siepv/siep/data/TS_All_Tissues/h5ad")


files <- list.files(pattern = "\\.h5ad$")
for (file in files) {
    dest_file <- file.path("../h5seurat", gsub("\\.h5ad$", ".h5seurat", basename(file)))
    SeuratDisk::Convert(file, dest = dest_file, overwrite = TRUE)
}
