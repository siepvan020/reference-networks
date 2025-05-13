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


# Function to merge cell types into a new label
merge_cell_types <- function(object, cell_type_mappings) {
    for (new_label in names(cell_type_mappings)) {
        old_labels <- cell_type_mappings[[new_label]]
        object@meta.data$cell_type[object@meta.data$cell_type %in% old_labels] <- new_label
    }
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


# Set active identity to cell type
Idents(blood_object_filter) <- blood_object_filter@meta.data$cell_type
Idents(lung_object_filter) <- lung_object_filter@meta.data$cell_type
Idents(fat_object_filter) <- fat_object_filter@meta.data$cell_type
Idents(kidney_object_filter) <- kidney_object_filter@meta.data$cell_type
Idents(liver_object_filter) <- liver_object_filter@meta.data$cell_type


#### 5. Convert ensembl ids to gene names and remove duplicates ####

# Extract feature names and Ensembl IDs
features <- data.frame(
    gene_name = gsub("_ENSG[0-9]+", "", as.character(blood_object_filter@assays$RNA@meta.features$feature_name)),
    ensembl_id = as.character(blood_object_filter@assays$RNA@counts@Dimnames[[1]])
)

# Identify duplicate gene names
duplicates <- features[duplicated(features$gene_name) | duplicated(features$gene_name, fromLast = TRUE), ]

# Get unique Ensembl IDs of duplicates
genes_to_exclude <- unique(duplicates$ensembl_id)

# Filter out the duplicate genes from the count matrix
keep_genes <- !blood_object_filter@assays$RNA@counts@Dimnames[[1]] %in% genes_to_exclude

# Update the gene names and filter the expression matrix accordingly
blood_object_filter@assays$RNA@counts <- blood_object_filter@assays$RNA@counts[keep_genes, ]
blood_object_filter@assays$RNA@data <- blood_object_filter@assays$RNA@counts[keep_genes, ]
blood_object_filter@assays$RNA@counts@Dimnames[[1]] <- features$gene_name[keep_genes]
blood_object_filter@assays$RNA@data@Dimnames[[1]] <- features$gene_name[keep_genes]


# Apply the same filtering to the lung object
lung_object_filter@assays$RNA@counts <- lung_object_filter@assays$RNA@counts[keep_genes, ]
lung_object_filter@assays$RNA@data <- lung_object_filter@assays$RNA@counts[keep_genes, ]
lung_object_filter@assays$RNA@counts@Dimnames[[1]] <- features$gene_name[keep_genes]
lung_object_filter@assays$RNA@data@Dimnames[[1]] <- features$gene_name[keep_genes]




setwd("/div/pythagoras/u1/siepv/siep/Analysis_v2/output/oplossing!!!!/prep")

# Save the filtered objects
saveRDS(blood_object_filter, "blood_filter.rds")
saveRDS(lung_object_filter, "lung_filter.rds")
saveRDS(fat_object_filter, "fat_filter.rds")
saveRDS(kidney_object_filter, "kidney_filter.rds")
saveRDS(liver_object_filter, "liver_filter.rds")
