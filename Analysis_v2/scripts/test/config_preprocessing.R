#!/usr/bin/env Rscript


#### 1. Setup ####

progress_file <- "/div/pythagoras/u1/siepv/siep/Analysis_v2/output/log/preprocess.log"

# Start fresh log or append
cat("- Starting script at -", format(Sys.time()), "\n", file = progress_file)

if (!requireNamespace("SeuratObject", quietly = TRUE)) remotes::install_version("SeuratObject", "4.1.4", repos = c("https://satijalab.r-universe.dev", getOption("repos")))
if (!requireNamespace("Seurat", quietly = TRUE)) remotes::install_version("Seurat", "4.4.0", repos = c("https://satijalab.r-universe.dev", getOption("repos")))

library(Seurat)
library(SeuratObject)
library(biomaRt)
library(dplyr)
library(ggplot2)
library(sva)

library(yaml)
library(glue)
library(purrr)

# Set working directory
setwd("/div/pythagoras/u1/siepv/siep/Analysis_v2")

# Load config file
config <- yaml::read_yaml("data/config/config.yaml")


#### 2. Load data ####

# Convert NULL to empty character vector
.null_to_chr <- function(x) if (is.null(x)) character(0) else x

# Reads one tissue object, merges and corrects celltype annotations, then
# keeps only the desired cell types and returns the filtered Seurat object.
process_tissue <- function(spec, tissue) {
    object <- readRDS(glue("data/rds/{spec$file}"))

    mappings <- map(spec$cells, .null_to_chr) |> # ensure vectors
        keep(~ length(.x) > 0) # drop empty ones

    object <- merge_cell_types(object, mappings) # Returns object with updated cell types in metadata
    object <- object[, object@meta.data$cell_type %in% names(spec$cells)]
    object
}

# Function to merge cell types into a new label
merge_cell_types <- function(object, cell_type_mappings) {
    object <- Seurat::DietSeurat(object)

    object@meta.data$cell_type <- as.character(object@meta.data$cell_type)

    for (new_label in names(cell_type_mappings)) {
        old_labels <- cell_type_mappings[[new_label]]
        object@meta.data$cell_type[object@meta.data$cell_type %in% old_labels] <- new_label
    }
    object@meta.data$cell_type <- as.factor(object@meta.data$cell_type)
    Idents(object) <- "cell_type"

    return(object)
}

objects <- imap(config, process_tissue)
