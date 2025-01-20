#!/usr/bin/env Rscript

# Set CRAN mirror
options(repos = c(CRAN = "https://cloud.r-project.org"))


# Load required packages
install.packages("remotes")
remotes::install_github("mojaveazure/seurat-disk")
install.packages("SCORPION")
install.packages("hdf5r", dependencies = TRUE)


library(SeuratDisk, verbose = TRUE)
library(SCORPION)
library(Seurat)
library(ggplot2)
