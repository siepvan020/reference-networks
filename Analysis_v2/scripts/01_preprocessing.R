#!/usr/bin/env Rscript

# Load required packages
if (!requireNamespace("Seurat", quietly = TRUE)) install.packages("Seurat")
# if (!requireNamespace("harmony", quietly = TRUE)) install.packages("harmony")

library(Seurat)
# library(harmony)

# Set working directory
setwd("/div/pythagoras/u1/siepv/siep/Analysis_v2/data/rds")


#### 1. Load data ####
blood_object <- readRDS("TS_v2_Blood.rds")
lung_object <- readRDS("TS_v2_Lung.rds")
print("Loaded data")


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

print("Filtered data")

#### 4. Batch correction ####




#### 5. Save preprocessed data ####
setwd("/div/pythagoras/u1/siepv/siep/Analysis_v2/output/preprocessing")

saveRDS(blood_object, file = "blood_preprocessed.rds")
saveRDS(lung_object, file = "lung_preprocessed.rds")

print("Saved filtered, batch corrected objects")
