#!/usr/bin/env Rscript


if (!requireNamespace("SeuratObject", quietly = TRUE)) remotes::install_version("SeuratObject", "4.1.4", repos = c("https://satijalab.r-universe.dev", getOption("repos")))
if (!requireNamespace("Seurat", quietly = TRUE)) remotes::install_version("Seurat", "4.4.0", repos = c("https://satijalab.r-universe.dev", getOption("repos")))

library(Seurat)
library(SeuratObject)

# Set working directory
setwd("/div/pythagoras/u1/siepv/siep/Analysis_v2/data")


#### 1. Load data ####
blood_object <- readRDS("rds/TS_v2_Blood.rds")
lung_object <- readRDS("rds/TS_v2_Lung.rds")
fat_object <- readRDS("rds/TS_v2_Fat.rds")
kidney_object <- readRDS("rds/TS_v2_Kidney.rds")
liver_object <- readRDS("rds/TS_v2_Liver.rds")

tf <- read.delim("priors/motif_prior_names_2024.tsv", header = FALSE, sep = "\t")


#### 2. Correct cell type annotations ####

# Function to merge cell types into a new label
merge_cell_types <- function(object, cell_type_mappings) {
    object <- Seurat::DietSeurat(object, dimreducs = c("pca", "umap"))

    object@meta.data$cell_type <- as.character(object@meta.data$cell_type)

    for (new_label in names(cell_type_mappings)) {
        old_labels <- cell_type_mappings[[new_label]]
        object@meta.data$cell_type[object@meta.data$cell_type %in% old_labels] <- new_label
    }
    object@meta.data$cell_type <- as.factor(object@meta.data$cell_type)
    Idents(object) <- "cell_type"

    return(object)
}

# Define cell type mappings for each object
blood_mappings <- list(
    "monocyte" = c("classical monocyte", "non-classical monocyte", "intermediate monocyte", "monocyte"),
    "cd4_positive_t_cell" = c("CD4-positive, alpha-beta T cell", "naive thymus-derived CD4-positive, alpha-beta T cell"),
    "cd8_positive_t_cell" = c("CD8-positive, alpha-beta T cell")
)

lung_mappings <- list(
    "monocyte" = c("classical monocyte", "non-classical monocyte", "intermediate monocyte", "monocyte"),
    "cd4_positive_t_cell" = c("CD4-positive, alpha-beta T cell"),
    "cd8_positive_t_cell" = c("CD8-positive, alpha-beta T cell"),
    "type_ii_pneumocyte" = c("pulmonary alveolar type 2 cell")
)

fat_mappings <- list(
    "monocyte" = c("classical monocyte", "non-classical monocyte", "intermediate monocyte", "monocyte"),
    "cd4_positive_t_cell" = c("CD4-positive, alpha-beta T cell"),
    "cd8_positive_t_cell" = c("CD8-positive, alpha-beta T cell"),
    "adipose_stem_cell" = c("mesenchymal stem cell of adipose tissue")
)

kidney_mappings <- list(
    "monocyte" = c("non-classical monocyte", "intermediate monocyte", "monocyte"),
    "cd4_positive_t_cell" = c("CD4-positive, alpha-beta T cell"),
    "cd8_positive_t_cell" = c("CD8-positive, alpha-beta T cell"),
    "kidney_epithelial_cell" = c("kidney epithelial cell")
)

liver_mappings <- list(
    "monocyte" = c("classical monocyte", "non-classical monocyte", "intermediate monocyte", "monocyte"),
    "cd4_positive_t_cell" = c("CD4-positive, alpha-beta T cell"),
    "cd8_positive_t_cell" = c("CD8-positive, alpha-beta T cell")
)

# Apply the function to each object
blood_object <- merge_cell_types(blood_object, blood_mappings)
lung_object <- merge_cell_types(lung_object, lung_mappings)
fat_object <- merge_cell_types(fat_object, fat_mappings)
kidney_object <- merge_cell_types(kidney_object, kidney_mappings)
liver_object <- merge_cell_types(liver_object, liver_mappings)


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


#### 4. Convert ensembl ids to gene names and remove duplicates ####

# Extract feature names and Ensembl IDs, filtering out gene names matching the regex
features <- data.frame(
    gene_name = gsub("_ENSG[0-9]+", "", as.character(blood_object_filter@assays$RNA@meta.features$feature_name)),
    ensembl_id = as.character(rownames(blood_object_filter))
)

# Identify duplicate gene names
duplicates <- features[duplicated(features$gene_name) | duplicated(features$gene_name, fromLast = TRUE), ]

# Get unique Ensembl IDs of duplicates
genes_to_exclude <- unique(duplicates$ensembl_id)

# Exclude duplicates
keep_genes <- !rownames(blood_object_filter) %in% genes_to_exclude

# Function to update gene names and filter the expression matrix
update_gene_names_and_filter <- function(object, keep_genes, features) {

    # Subset genes that are not mitochondrial
    non_mt_genes <- !object@assays$RNA@meta.features$mt
    keep_genes <- non_mt_genes & keep_genes

    # Update gene names and filter duplicates
    filtered_object <- object[keep_genes, ]

    # List of new gene names
    new_names <- gsub("\\.[0-9]+$", "", features$gene_name[keep_genes])

    assayobj <- filtered_object@assays$RNA
    rownames(assayobj@meta.features) <- new_names

    matrix_n <- SeuratObject::GetAssayData(assayobj)

    matrix_n@Dimnames[[1]] <- new_names

    assayobj@counts <- matrix_n
    assayobj@data <- matrix_n

    object@assays[["RNA"]] <- assayobj

    return(object)
}

# Apply the function to each object
blood_object_filter <- update_gene_names_and_filter(blood_object_filter, keep_genes, features)
lung_object_filter <- update_gene_names_and_filter(lung_object_filter, keep_genes, features)
fat_object_filter <- update_gene_names_and_filter(fat_object_filter, keep_genes, features)
kidney_object_filter <- update_gene_names_and_filter(kidney_object_filter, keep_genes, features)
liver_object_filter <- update_gene_names_and_filter(liver_object_filter, keep_genes, features)

setwd("/div/pythagoras/u1/siepv/siep/Analysis_v2/output/preprocessing")

# Save the filtered objects
saveRDS(blood_object_filter, "blood_filter.rds")
saveRDS(lung_object_filter, "lung_filter.rds")
saveRDS(fat_object_filter, "fat_filter.rds")
saveRDS(kidney_object_filter, "kidney_filter.rds")
saveRDS(liver_object_filter, "liver_filter.rds")
