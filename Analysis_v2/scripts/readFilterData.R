#!/usr/bin/env Rscript

# Load required packages
if (!requireNamespace("Seurat", quietly = TRUE)) install.packages("Seurat")
if (!requireNamespace("harmony", quietly = TRUE)) install.packages("harmony")

library(Seurat)
library(harmony)

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
Idents(filtered_blood) <- filtered_blood@meta.data$cell_type

table(filtered_blood@meta.data$cell_type)


setwd("/div/pythagoras/u1/siepv/siep/Analysis_v2/output/preprocessing")
lung_filter_obj <- readRDS("lung_preprocessed.rds")
lung_batch_obj <- readRDS("lung_preprocessed_batch.rds")

lung_filter_obj@reductions <- list()
lung_batch_obj@reductions <- list()

lung_filter_obj <- FindVariableFeatures(lung_filter_obj)
lung_filter_obj <- ScaleData(lung_filter_obj)
lung_filter_obj <- RunPCA(lung_filter_obj, verbose = FALSE)
lung_filter_obj <- RunUMAP(lung_filter_obj, dims = 1:30)

lung_batch_obj <- FindVariableFeatures(lung_batch_obj)
lung_batch_obj <- ScaleData(lung_batch_obj)
lung_batch_obj <- RunPCA(lung_batch_obj, verbose = FALSE)
lung_batch_obj <- RunUMAP(lung_batch_obj, dims = 1:30)

p1 <- DimPlot(lung_filter_obj, reduction = "umap", group.by = "cell_type")
p2 <- DimPlot(lung_batch_obj, reduction = "umap", group.by = "cell_type")
p1 + p2

p3 <- DimPlot(lung_filter_obj, reduction = "umap", group.by = "donor_id")
p4 <- DimPlot(lung_batch_obj, reduction = "umap", group.by = "donor_id")
p3 + p4


batch_blood_obj <- harmony::RunHarmony(filtered_blood, group.by.vars = "donor_id")
DimPlot(batch_blood_obj, reduction = "harmony", group.by = "cell_type")
DimPlot(batch_blood_obj, reduction = "harmony", group.by = "donor_id")


setwd("/div/pythagoras/u1/siepv/siep/Analysis_v2/output/preprocessing")

saveRDS(batch_blood, file = "blood_preprocessed.rds")
print("Saved filtered, batch corrected blood object")





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
Idents(filtered_lung) <- filtered_lung@meta.data$cell_type

table(filtered_lung@meta.data$cell_type)

# batch_lung <- harmony::RunHarmony(filtered_lung, group.by.vars = "donor_id")

setwd("/div/pythagoras/u1/siepv/siep/Analysis_v2/output/preprocessing")

saveRDS(batch_lung, file = "lung_preprocessed.rds")
print("Saved filtered, batch corrected lung object")
