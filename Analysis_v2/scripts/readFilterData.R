#!/usr/bin/env Rscript

# Load required packages
if (!requireNamespace("Seurat", quietly = TRUE)) install.packages("Seurat")

library(Seurat)

# Set working directory
setwd("/div/pythagoras/u1/siepv/siep/Analysis_v2/data/rds")

# ------------------------------ Blood ------------------------------

blood_object <- readRDS("TS_v2_Blood.rds")

# Convert cell type to character
blood_object@meta.data$cell_type <- as.character(blood_object@meta.data$cell_type)

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

# Define which cell types should be kept
blood_types_to_keep <- c(
    "cd4_positive_t_cell",
    "cd8_positive_t_cell",
    "monocyte",
    "platelet"
)

# Get the barcodes (rownames) of the cells that match the selected types
blood_to_keep <- rownames(blood_object@meta.data[blood_object@meta.data$cell_type %in% blood_types_to_keep, ])

# Filter the Seurat object based on selected cell barcodes
filtered_blood <- blood_object[, blood_to_keep]

table(filtered_blood@meta.data$cell_type)

setwd("/div/pythagoras/u1/siepv/siep/Analysis_v2/output/celltype_filtering")

saveRDS(filtered_blood, file = "blood_celltypes.rds")
print("Saved filtered blood cell types")

# ------------------------------ Lung ------------------------------

setwd("/div/pythagoras/u1/siepv/siep/Analysis_v2/data/rds")

lung_object <- readRDS("TS_v2_Lung.rds")

# Convert cell type to character
lung_object@meta.data$cell_type <- as.character(lung_object@meta.data$cell_type)

# Define which cell types should be merged into a new label
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


# Define which cell types should be kept
lung_types_to_keep <- c(
    "cd4_positive_t_cell",
    "cd8_positive_t_cell",
    "monocyte",
    "type_ii_pneumocyte"
)

# Get the barcodes (rownames) of the cells that match the selected types
lung_to_keep <- rownames(lung_object@meta.data[lung_object@meta.data$cell_type %in% lung_types_to_keep, ])

# Filter the Seurat object based on selected cell barcodes
filtered_lung <- lung_object[, lung_to_keep]

table(filtered_lung@meta.data$cell_type)

setwd("/div/pythagoras/u1/siepv/siep/Analysis_v2/output/celltype_filtering")

saveRDS(filtered_lung, file = "lung_celltypes.rds")
print("Saved filtered lung cell types")
