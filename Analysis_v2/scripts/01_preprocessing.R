#!/usr/bin/env Rscript

library(Seurat)

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
    # blood_object@meta.data$cell_type <- as.character(blood_object@meta.data$cell_type)
    # object@assays$RNA@data <- NULL
    # for (new_label in names(cell_type_mappings)) {
    #     old_labels <- cell_type_mappings[[new_label]]
    #     object@meta.data$cell_type[object@meta.data$cell_type %in% old_labels] <- new_label
    # }

    # object <- Seurat::DietSeurat(object, layers="counts", dimreducs=c("pca", "umap"))

    idents.df = data.frame("cell_type" = object$cell_type, "donor_id" = object$donor_id)
    export = CreateSeuratObject(counts = LayerData(object, assay = "RNA", layer = "counts"), meta.data = idents.df)


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
    keep_genes <- keep_genes & non_mt_genes

    # Update gene names and filter duplicates
    filtered_counts <- object[["RNA"]]$counts[keep_genes, ]
    rownames(filtered_counts) <- features$gene_name[keep_genes]
    rownames(filtered_counts) <- gsub("\\.[0-9]+$", "", rownames(filtered_counts))

    # Filter out rows with ENSG gene names
    # filtered_counts <- filtered_counts[!grepl("^ENSG[0-9]+(\\.[0-9]+)?$", rownames(filtered_counts)), ]

    # Only take the intersect between priors and gex genes
    common_genes <- intersect(rownames(filtered_counts), unique(tf$V2))
    filtered_counts <- filtered_counts[common_genes, ]

    # Update the object's counts
    object[["RNA"]]$counts <- filtered_counts
    
    # object[["RNA"]] <- Seurat::CreateAssayObject(counts = object[["RNA"]]@counts)

    # object[["RNA"]]$data <- Seurat::NormalizeData(object[["RNA"]]$counts)

    # Ensure the data layer contains the same genes as the counts layer
    # filtered_data <- object@assays$RNA@data[rownames(filtered_counts), , drop = FALSE]
    # object[["RNA"]]$data <- NULL#filtered_data

    return(object)
}

#
#
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




#### 5. Batch correction ####
# note that kidney is not included, as this tissue only has one donor

# sce_blood <- Seurat::as.SingleCellExperiment(blood_object_filter)
# sce_lung <- Seurat::as.SingleCellExperiment(lung_object_filter)
# sce_fat <- Seurat::as.SingleCellExperiment(fat_object_filter)
# sce_liver <- Seurat::as.SingleCellExperiment(liver_object_filter)

# cat("Run ComBatSeq on blood", format(Sys.time()), "\n", file = progress_file, append = TRUE)
# sce_blood <- singleCellTK::runComBatSeq(sce_blood, useAssay = "counts", batch = "donor_id", assayName = "ComBatSeq")

# cat("Run ComBatSeq on lung", format(Sys.time()), "\n\n", file = progress_file, append = TRUE)
# sce_lung <- singleCellTK::runComBatSeq(sce_lung, useAssay = "counts", batch = "donor_id", assayName = "ComBatSeq")

# cat("Run ComBatSeq on fat", format(Sys.time()), "\n\n", file = progress_file, append = TRUE)
# sce_fat <- singleCellTK::runComBatSeq(sce_fat, useAssay = "counts", batch = "donor_id", assayName = "ComBatSeq")

# cat("Run ComBatSeq on liver", format(Sys.time()), "\n\n", file = progress_file, append = TRUE)
# sce_liver <- singleCellTK::runComBatSeq(sce_liver, useAssay = "counts", batch = "donor_id", assayName = "ComBatSeq")

# cat("Assays in SCE Blood:", assayNames(sce_blood), "\n", file = progress_file, append = TRUE)
# cat("Assays in SCE Lung:", assayNames(sce_lung), "\n\n", file = progress_file, append = TRUE)
# cat("Assays in SCE Fat:", assayNames(sce_fat), "\n\n", file = progress_file, append = TRUE)
# cat("Assays in SCE Liver:", assayNames(sce_liver), "\n\n", file = progress_file, append = TRUE)

# cat("Save SCE objects at", format(Sys.time()), "\n\n\n", file = progress_file, append = TRUE)
# saveRDS(sce_blood, file = "blood_sce.rds")
# saveRDS(sce_lung, file = "lung_sce.rds")
# saveRDS(sce_fat, file = "fat_sce.rds")
# saveRDS(sce_liver, file = "liver_sce.rds")

# cat("Convert objects back to Seurat at", format(Sys.time()), "\n", file = progress_file, append = TRUE)
# batch_blood_obj <- as.Seurat(sce_blood)
# batch_lung_obj <- as.Seurat(sce_lung)
# batch_fat_obj <- as.Seurat(sce_fat)
# batch_liver_obj <- as.Seurat(sce_liver)

# cat("Save batch corrected seurat object at", format(Sys.time()), "\n\n\n", file = progress_file, append = TRUE)
# saveRDS(batch_blood_obj, file = "blood_filter_batch.rds")
# saveRDS(batch_lung_obj, file = "lung_filter_batch.rds")
# saveRDS(batch_fat_obj, file = "fat_filter_batch.rds")
# saveRDS(batch_liver_obj, file = "liver_filter_batch.rds")

# # Ensure the corrected matrix is assigned to RNA counts in Seurat
# update_counts_with_combatseq <- function(sce_object, seurat_object, tissue_name) {
#     if ("ComBatSeq" %in% assayNames(sce_object)) {
#         seurat_object@assays$RNA@counts <- assay(sce_object, "ComBatSeq")
#     } else {
#         cat(paste("Warning: ComBatSeq assay missing in SCE", tissue_name, "object!\n"), file = progress_file, append = TRUE)
#     }
#     return(seurat_object)
# }

# batch_blood_obj <- update_counts_with_combatseq(sce_blood, batch_blood_obj, "blood")
# batch_lung_obj <- update_counts_with_combatseq(sce_lung, batch_lung_obj, "lung")
# batch_fat_obj <- update_counts_with_combatseq(sce_fat, batch_fat_obj, "fat")
# batch_liver_obj <- update_counts_with_combatseq(sce_liver, batch_liver_obj, "liver")


# #### 6. Save preprocessed data ####
# setwd("/div/pythagoras/u1/siepv/siep/Analysis_v2/output/preprocessing")

# cat("Overwrite batch corrected seurat object at", format(Sys.time()), "\n\n", file = progress_file, append = TRUE)
# saveRDS(batch_blood_obj, file = "blood_filter_batch.rds")
# saveRDS(batch_lung_obj, file = "lung_filter_batch.rds")
# saveRDS(batch_fat_obj, file = "fat_filter_batch.rds")
# saveRDS(batch_liver_obj, file = "liver_filter_batch.rds")

# cat("Saved filtered, batch corrected objects at", format(Sys.time()), "\n", file = progress_file, append = TRUE)
