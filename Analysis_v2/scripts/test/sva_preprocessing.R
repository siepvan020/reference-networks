#!/usr/bin/env Rscript

progress_file <- "/div/pythagoras/u1/siepv/siep/Analysis_v2/output/log/sva_preprocess.log"
cat("Starting script at", format(Sys.time()), "\n", file = progress_file)

# Load required packages
if (!requireNamespace("devtools", quietly = TRUE)) install.packages("devtools")

if (!requireNamespace("sva", quietly = TRUE)) devtools::install_github("zhangyuqing/sva-devel")
if (!requireNamespace("Seurat", quietly = TRUE)) install.packages("Seurat")

if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
# if (!requireNamespace("singleCellTK", quietly = TRUE)) BiocManager::install("singleCellTK")
if (!requireNamespace("Matrix", quietly = TRUE)) install.packages("Matrix")
if (!requireNamespace("ggplot2", quietly = TRUE)) install.packages("ggplot2")

library(Seurat)
library(sva)
library(Matrix)
# library(singleCellTK)
library(ggplot2)
# library(harmony)

cat("Loaded required packages at", format(Sys.time()), "\n", file = progress_file, append = TRUE)


# Set working directory
setwd("/div/pythagoras/u1/siepv/siep/Analysis_v2/data/rds")


#### 1. Load data ####
blood_object <- readRDS("TS_v2_Blood.rds")
lung_object <- readRDS("TS_v2_Lung.rds")
cat("Loaded data at", format(Sys.time()), "\n", file = progress_file, append = TRUE)


#### 2. Correct cell type annotations ####
blood_object@meta.data$cell_type <- as.character(blood_object@meta.data$cell_type)
lung_object@meta.data$cell_type <- as.character(lung_object@meta.data$cell_type)

# Define which cell types should be merged into a new label
blood_object@meta.data$cell_type[blood_object@meta.data$cell_type %in% c(
    "classical monocyte",
    "non-classical monocyte",
    "intermediate monocyte",
    "monocyte"
)] <- "monocyte"
blood_object@meta.data$cell_type[blood_object@meta.data$cell_type %in% c(
    "CD4-positive, alpha-beta T cell",
    "naive thymus-derived CD4-positive, alpha-beta T cell"
)] <- "cd4_positive_t_cell"
blood_object@meta.data$cell_type[blood_object@meta.data$cell_type %in% c(
    "CD8-positive, alpha-beta T cell"
)] <- "cd8_positive_t_cell"

lung_object@meta.data$cell_type[lung_object@meta.data$cell_type %in% c(
    "classical monocyte",
    "non-classical monocyte",
    "intermediate monocyte",
    "monocyte"
)] <- "monocyte"
lung_object@meta.data$cell_type[lung_object@meta.data$cell_type %in% c(
    "CD4-positive, alpha-beta T cell"
)] <- "cd4_positive_t_cell"
lung_object@meta.data$cell_type[lung_object@meta.data$cell_type %in% c(
    "CD8-positive, alpha-beta T cell"
)] <- "cd8_positive_t_cell"
lung_object@meta.data$cell_type[lung_object@meta.data$cell_type %in% c(
    "pulmonary alveolar type 2 cell"
)] <- "type_ii_pneumocyte"


#### 3. Filter only relevant cell types ####
cell_types_to_keep <- c("cd4_positive_t_cell", "cd8_positive_t_cell", "monocyte", "platelet", "type_ii_pneumocyte")
cells_to_keep_blood <- rownames(blood_object@meta.data[blood_object@meta.data$cell_type %in% cell_types_to_keep, ])
cells_to_keep_lung <- rownames(lung_object@meta.data[lung_object@meta.data$cell_type %in% cell_types_to_keep, ])

# Apply filtering
blood_object <- blood_object[, cells_to_keep_blood]
lung_object <- lung_object[, cells_to_keep_lung]
Idents(blood_object) <- blood_object@meta.data$cell_type
Idents(lung_object) <- lung_object@meta.data$cell_type
cat("Filtered data at", format(Sys.time()), "\n", file = progress_file, append = TRUE)

setwd("/div/pythagoras/u1/siepv/siep/Analysis_v2/output/sva_preprocessing")
saveRDS(blood_object, file = "blood_preprocessed.rds")
saveRDS(lung_object, file = "lung_preprocessed.rds")
cat("Saved filtered data at", format(Sys.time()), "\n\n\n", file = progress_file, append = TRUE)


#### 4. Batch correction ####

cat("Run ComBatSeq on blood", format(Sys.time()), "\n", file = progress_file, append = TRUE)
blood_input_matrix <- blood_object@assays$RNA@counts
batch_matrix_blood <- sva::ComBat_seq(counts = blood_input_matrix, batch = blood_object@meta.data$donor_id)
batch_blood_obj <- blood_object
batch_blood_obj@assays$RNA@counts <- as(batch_matrix_blood, "dgCMatrix")
saveRDS(batch_blood_obj, file = "blood_preprocessed_batch.rds")


cat("Run ComBatSeq on lung", format(Sys.time()), "\n", file = progress_file, append = TRUE)
lung_input_matrix <- lung_object@assays$RNA@counts
batch_matrix_lung <- sva::ComBat_seq(counts = lung_input_matrix, batch = lung_object@meta.data$donor_id)
batch_lung_obj <- lung_object
batch_lung_obj@assays$RNA@counts <- as(batch_matrix_lung, "dgCMatrix")
saveRDS(batch_lung_obj, file = "lung_preprocessed_batch.rds")


#### 5. Save preprocessed data ####
setwd("/div/pythagoras/u1/siepv/siep/Analysis_v2/output/preprocessing")


cat("Save preprocessed data at", format(Sys.time()), "\n", file = progress_file, append = TRUE)
saveRDS(blood_object, file = "blood_preprocessed.rds")
saveRDS(lung_object, file = "lung_preprocessed.rds")


cat("Saved filtered, batch corrected objects at", format(Sys.time()), "\n", file = progress_file, append = TRUE)
