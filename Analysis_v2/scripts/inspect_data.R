#!/usr/bin/env Rscript

# Load required packages
if (!requireNamespace("Seurat", quietly = TRUE)) install.packages("Seurat")

library(Seurat)

# Set working directory
setwd("/div/pythagoras/u1/siepv/siep/Analysis_v2/data/rds")

# Load objects
blood_object <- readRDS("TS_v2_Blood.rds")
lung_object <- readRDS("TS_v2_Lung.rds")
