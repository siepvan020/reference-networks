#!/usr/bin/env Rscript

library(Seurat)


# Set working directory
setwd("/div/pythagoras/u1/siepv/siep/Analysis_v2/data/rds")


#### 1. Load data ####
blood_object <- readRDS("TS_v2_Blood.rds")
lung_object <- readRDS("TS_v2_Lung.rds")
fat_object <- readRDS("TS_v2_Fat.rds")
kidney_object <- readRDS("TS_v2_Kidney.rds")
liver_object <- readRDS("TS_v2_Liver.rds")


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
blood_object_filter <- blood_object[, cells_to_keep_blood]
lung_object_filter <- lung_object[, cells_to_keep_lung]
fat_object_filter <- fat_object[, cells_to_keep_fat]
kidney_object_filter <- kidney_object[, cells_to_keep_kidney]
liver_object_filter <- liver_object[, cells_to_keep_liver]


# Not really needed
Idents(blood_object_filter) <- blood_object_filter@meta.data$cell_type
Idents(lung_object_filter) <- lung_object_filter@meta.data$cell_type
Idents(fat_object_filter) <- fat_object_filter@meta.data$cell_type
Idents(kidney_object_filter) <- kidney_object_filter@meta.data$cell_type
Idents(liver_object_filter) <- liver_object_filter@meta.data$cell_type


# Check how many genes are not expressed in any cell
table(rowSums(blood_object_filter@assays$RNA@counts > 0) == 0)
table(rowSums(lung_object_filter@assays$RNA@counts > 0) == 0)
table(rowSums(fat_object_filter@assays$RNA@counts > 0) == 0)
table(rowSums(kidney_object_filter@assays$RNA@counts > 0) == 0)
table(rowSums(liver_object_filter@assays$RNA@counts > 0) == 0)
