#!/usr/bin/env Rscript


progress_file <- "/div/pythagoras/u1/siepv/siep/Analysis_v2/output/log/preprocess.log"
cat("Starting script at", format(Sys.time()), "\n", file = progress_file)

# Load required packages
if (!requireNamespace("devtools", quietly = TRUE)) install.packages("devtools")
if (!requireNamespace("Seurat", quietly = TRUE)) install.packages("Seurat")
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
if (!requireNamespace("singleCellTK", quietly = TRUE)) BiocManager::install("singleCellTK")
if (!requireNamespace("Matrix", quietly = TRUE)) install.packages("Matrix")
if (!requireNamespace("ggplot2", quietly = TRUE)) install.packages("ggplot2")

library(Seurat)
library(Matrix)
library(singleCellTK)
library(ggplot2)

cat("Loaded required packages at", format(Sys.time()), "\n", file = progress_file, append = TRUE)



# Set working directory
setwd("/div/pythagoras/u1/siepv/siep/Analysis_v2/data/rds")


#### 1. Load data ####
blood_object <- readRDS("TS_v2_Blood.rds")
lung_object <- readRDS("TS_v2_Lung.rds")
fat_object <- readRDS("TS_v2_Fat.rds")
kidney_object <- readRDS("TS_v2_Kidney.rds")
liver_object <- readRDS("TS_v2_Liver.rds")
cat("Loaded data at", format(Sys.time()), "\n", file = progress_file, append = TRUE)


#### 2. Correct cell type annotations ####
blood_object@meta.data$cell_type <- as.character(blood_object@meta.data$cell_type)
lung_object@meta.data$cell_type <- as.character(lung_object@meta.data$cell_type)
fat_object@meta.data$cell_type <- as.character(fat_object@meta.data$cell_type)
kidney_object@meta.data$cell_type <- as.character(kidney_object@meta.data$cell_type)
liver_object@meta.data$cell_type <- as.character(liver_object@meta.data$cell_type)


# Define which cell types should be merged into a new label
#################################################################################### Blood
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
#################################################################################### Lung
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
#################################################################################### Fat
fat_object@meta.data$cell_type[fat_object@meta.data$cell_type %in% c(
    "classical monocyte",
    "non-classical monocyte",
    "intermediate monocyte",
    "monocyte"
)] <- "monocyte"
fat_object@meta.data$cell_type[fat_object@meta.data$cell_type %in% c(
    "CD4-positive, alpha-beta T cell"
)] <- "cd4_positive_t_cell"
fat_object@meta.data$cell_type[fat_object@meta.data$cell_type %in% c(
    "CD8-positive, alpha-beta T cell"
)] <- "cd8_positive_t_cell"
fat_object@meta.data$cell_type[fat_object@meta.data$cell_type %in% c(
    "mesenchymal stem cell of adipose tissue"
)] <- "adipose_stem_cell"
#################################################################################### Kidney
kidney_object@meta.data$cell_type[kidney_object@meta.data$cell_type %in% c(
    "non-classical monocyte",
    "intermediate monocyte",
    "monocyte"
)] <- "monocyte"
kidney_object@meta.data$cell_type[kidney_object@meta.data$cell_type %in% c(
    "CD4-positive, alpha-beta T cell"
)] <- "cd4_positive_t_cell"
kidney_object@meta.data$cell_type[kidney_object@meta.data$cell_type %in% c(
    "CD8-positive, alpha-beta T cell"
)] <- "cd8_positive_t_cell"
kidney_object@meta.data$cell_type[kidney_object@meta.data$cell_type %in% c(
    "kidney epithelial cell"
)] <- "kidney_epithelial_cell"
#################################################################################### Liver
liver_object@meta.data$cell_type[liver_object@meta.data$cell_type %in% c(
    "classical monocyte",
    "non-classical monocyte",
    "intermediate monocyte",
    "monocyte"
)] <- "monocyte"
liver_object@meta.data$cell_type[liver_object@meta.data$cell_type %in% c(
    "CD4-positive, alpha-beta T cell"
)] <- "cd4_positive_t_cell"
liver_object@meta.data$cell_type[liver_object@meta.data$cell_type %in% c(
    "CD8-positive, alpha-beta T cell"
)] <- "cd8_positive_t_cell"


#### 3. Filter only relevant cell types ####
common_celltypes <- c("cd4_positive_t_cell", "cd8_positive_t_cell", "monocyte")
cells_to_keep_blood <- rownames(blood_object@meta.data[blood_object@meta.data$cell_type %in% append(common_celltypes, "platelet"), ])
cells_to_keep_lung <- rownames(lung_object@meta.data[lung_object@meta.data$cell_type %in% append(common_celltypes, "type_ii_pneumocyte"), ])
cells_to_keep_fat <- rownames(fat_object@meta.data[fat_object@meta.data$cell_type %in% append(common_celltypes, "adipose_stem_cell"), ])
cells_to_keep_kidney <- rownames(kidney_object@meta.data[kidney_object@meta.data$cell_type %in% append(common_celltypes, "kidney_epithelial_cell"), ])
cells_to_keep_liver <- rownames(liver_object@meta.data[liver_object@meta.data$cell_type %in% append(common_celltypes, "hepatocyte"), ])

# Apply filtering
blood_object <- blood_object[, cells_to_keep_blood]
lung_object <- lung_object[, cells_to_keep_lung]
fat_object <- fat_object[, cells_to_keep_fat]
kidney_object <- kidney_object[, cells_to_keep_kidney]
liver_object <- liver_object[, cells_to_keep_liver]

Idents(blood_object) <- blood_object@meta.data$cell_type
Idents(lung_object) <- lung_object@meta.data$cell_type
Idents(fat_object) <- fat_object@meta.data$cell_type
Idents(kidney_object) <- kidney_object@meta.data$cell_type
Idents(liver_object) <- liver_object@meta.data$cell_type
cat("Filtered data at", format(Sys.time()), "\n", file = progress_file, append = TRUE)

setwd("/div/pythagoras/u1/siepv/siep/Analysis_v2/output/preprocessing")
saveRDS(blood_object, file = "blood_filter.rds")
saveRDS(lung_object, file = "lung_filter.rds")
saveRDS(fat_object, file = "fat_filter.rds")
saveRDS(kidney_object, file = "kidney_filter.rds")
saveRDS(liver_object, file = "liver_filter.rds")
cat("Saved filtered data at", format(Sys.time()), "\n\n\n", file = progress_file, append = TRUE)


#### 4. Batch correction ####   note that kidney is not included, as this tissue only has one donor

sce_blood <- Seurat::as.SingleCellExperiment(blood_object)
sce_lung <- Seurat::as.SingleCellExperiment(lung_object)
sce_fat <- Seurat::as.SingleCellExperiment(fat_object)
sce_liver <- Seurat::as.SingleCellExperiment(liver_object)

cat("Run ComBatSeq on blood", format(Sys.time()), "\n", file = progress_file, append = TRUE)
sce_blood <- singleCellTK::runComBatSeq(sce_blood, useAssay = "counts", batch = "donor_id", assayName = "ComBatSeq")

cat("Run ComBatSeq on lung", format(Sys.time()), "\n\n", file = progress_file, append = TRUE)
sce_lung <- singleCellTK::runComBatSeq(sce_lung, useAssay = "counts", batch = "donor_id", assayName = "ComBatSeq")

cat("Run ComBatSeq on fat", format(Sys.time()), "\n\n", file = progress_file, append = TRUE)
sce_fat <- singleCellTK::runComBatSeq(sce_fat, useAssay = "counts", batch = "donor_id", assayName = "ComBatSeq")

cat("Run ComBatSeq on liver", format(Sys.time()), "\n\n", file = progress_file, append = TRUE)
sce_liver <- singleCellTK::runComBatSeq(sce_liver, useAssay = "counts", batch = "donor_id", assayName = "ComBatSeq")

cat("Assays in SCE Blood:", assayNames(sce_blood), "\n", file = progress_file, append = TRUE)
cat("Assays in SCE Lung:", assayNames(sce_lung), "\n\n", file = progress_file, append = TRUE)
cat("Assays in SCE Fat:", assayNames(sce_fat), "\n\n", file = progress_file, append = TRUE)
cat("Assays in SCE Liver:", assayNames(sce_liver), "\n\n", file = progress_file, append = TRUE)

cat("Save SCE objects at", format(Sys.time()), "\n\n\n", file = progress_file, append = TRUE)
saveRDS(sce_blood, file = "blood_sce.rds")
saveRDS(sce_lung, file = "lung_sce.rds")
saveRDS(sce_fat, file = "fat_sce.rds")
saveRDS(sce_liver, file = "liver_sce.rds")

cat("Convert objects back to Seurat at", format(Sys.time()), "\n", file = progress_file, append = TRUE)
batch_blood_obj <- as.Seurat(sce_blood)
batch_lung_obj <- as.Seurat(sce_lung)
batch_fat_obj <- as.Seurat(sce_fat)
batch_liver_obj <- as.Seurat(sce_liver)

cat("Save batch corrected seurat object at", format(Sys.time()), "\n\n\n", file = progress_file, append = TRUE)
saveRDS(batch_blood_obj, file = "blood_filter_batch.rds")
saveRDS(batch_lung_obj, file = "lung_filter_batch.rds")
saveRDS(batch_fat_obj, file = "fat_filter_batch.rds")
saveRDS(batch_liver_obj, file = "liver_filter_batch.rds")


# Ensure the corrected matrix is assigned to RNA counts in Seurat          ### CHANGE TO LESS REDUNDANCY
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
if ("ComBatSeq" %in% assayNames(sce_fat)) {
    batch_fat_obj@assays$RNA@counts <- assay(sce_fat, "ComBatSeq")
} else {
    cat("Warning: ComBatSeq assay missing in SCE fat object!\n", file = progress_file, append = TRUE)
}
if ("ComBatSeq" %in% assayNames(sce_liver)) {
    batch_liver_obj@assays$RNA@counts <- assay(sce_liver, "ComBatSeq")
} else {
    cat("Warning: ComBatSeq assay missing in SCE liver object!\n", file = progress_file, append = TRUE)
}


#### 5. Save preprocessed data ####
setwd("/div/pythagoras/u1/siepv/siep/Analysis_v2/output/preprocessing")


cat("Overwrite batch corrected seurat object at", format(Sys.time()), "\n\n", file = progress_file, append = TRUE)
saveRDS(batch_blood_obj, file = "blood_filter_batch.rds")
saveRDS(batch_lung_obj, file = "lung_filter_batch.rds")
saveRDS(batch_fat_obj, file = "fat_filter_batch.rds")
saveRDS(batch_liver_obj, file = "liver_filter_batch.rds")

cat("Saved filtered, batch corrected objects at", format(Sys.time()), "\n", file = progress_file, append = TRUE)
