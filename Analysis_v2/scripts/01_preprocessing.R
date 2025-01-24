#!/usr/bin/env Rscript



# Load required packages
if (!requireNamespace("devtools", quietly = TRUE)) install.packages("devtools")
if (!requireNamespace("Seurat", quietly = TRUE)) install.packages("Seurat")

if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
if (!requireNamespace("singleCellTK", quietly = TRUE)) BiocManager::install("singleCellTK")

# if (!requireNamespace("sva", quietly = TRUE)) devtools::install_github("zhangyuqing/sva-devel")
if (!requireNamespace("Matrix", quietly = TRUE)) install.packages("Matrix")
if (!requireNamespace("ggplot2", quietly = TRUE)) install.packages("ggplot2")

library(Seurat)
# library(sva)
library(Matrix)
library(singleCellTK)
library(ggplot2)
# library(harmony)

progress_file <- "/div/pythagoras/u1/siepv/siep/Analysis_v2/output/log/preprocess.log"

cat("Starting script at", format(Sys.time()), "\n", file = progress_file)


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

#### 4. Batch correction ####

sce_blood <- Seurat::as.SingleCellExperiment(blood_object)
sce_lung <- Seurat::as.SingleCellExperiment(lung_object)

cat("Run ComBatSeq on blood", format(Sys.time()), "\n", file = progress_file, append = TRUE)
sce_blood <- singleCellTK::runComBatSeq(sce_blood, useAssay = "counts", batch = "donor_id", assayName = "ComBatSeq")

cat("Run ComBatSeq on lung", format(Sys.time()), "\n", file = progress_file, append = TRUE)
sce_lung <- singleCellTK::runComBatSeq(sce_lung, useAssay = "counts", batch = "donor_id", assayName = "ComBatSeq")

cat("Convert objects back to Seurat at", format(Sys.time()), "\n", file = progress_file, append = TRUE)
batch_blood_obj <- as.Seurat(sce_blood)
batch_lung_obj <- as.Seurat(sce_lung)

# Ensure the corrected matrix is assigned to RNA counts in Seurat
if ("ComBatSeq" %in% assayNames(sce_blood)) {
    batch_blood_obj@assays$RNA@counts <- assay(sce_blood, "ComBatSeq")
} else {
    cat("Warning: ComBatSeq assay missing in SCE blood object!\n", file = progress_file, append = TRUE)
}

if ("ComBatSeq" %in% assayNames(sce_lung)) {
    batch_lung_obj@assays$RNA@counts <- assay(sce_lung, "ComBatSeq")
} else {
    cat("Warning: ComBatSeq assay missing in SCE lung object!\n", file = progress_file, append = TRUE)
}


#### 5. Save preprocessed data ####
setwd("/div/pythagoras/u1/siepv/siep/Analysis_v2/output/preprocessing")

cat("Save SCE objects at", format(Sys.time()), "\n", file = progress_file, append = TRUE)
saveRDS(sce_blood, file = "blood_sce.rds")
saveRDS(sce_lung, file = "lung_sce.rds")

cat("Save preprocessed data at", format(Sys.time()), "\n", file = progress_file, append = TRUE)
saveRDS(blood_object, file = "blood_preprocessed.rds")
saveRDS(lung_object, file = "lung_preprocessed.rds")

cat("Save batch corrected data at", format(Sys.time()), "\n", file = progress_file, append = TRUE)
saveRDS(batch_blood_obj, file = "blood_preprocessed_batch.rds")
saveRDS(batch_lung_obj, file = "lung_preprocessed_batch.rds")

cat("Saved filtered, batch corrected objects at", format(Sys.time()), "\n", file = progress_file, append = TRUE)
