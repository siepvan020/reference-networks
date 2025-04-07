#!/usr/bin/env Rscript


if (!requireNamespace("SeuratObject", quietly = TRUE)) remotes::install_version("SeuratObject", "4.1.4", repos = c("https://satijalab.r-universe.dev", getOption("repos")))
if (!requireNamespace("Seurat", quietly = TRUE)) remotes::install_version("Seurat", "4.4.0", repos = c("https://satijalab.r-universe.dev", getOption("repos")))

library(Seurat)
library(SeuratObject)
library(biomaRt)
library(dplyr)
library(ggplot2)
library(harmony)

# Set working directory
setwd("/div/pythagoras/u1/siepv/siep/Analysis_v2/data")


#### 1. Load data ####
blood_object <- readRDS("rds/TS_v2_Blood.rds")
lung_object <- readRDS("rds/TS_v2_Lung.rds")
fat_object <- readRDS("rds/TS_v2_Fat.rds")
kidney_object <- readRDS("rds/TS_v2_Kidney.rds")
liver_object <- readRDS("rds/TS_v2_Liver.rds")


#### 2. Correct cell type annotations ####

# Function to merge cell types into a new label
merge_cell_types <- function(object, cell_type_mappings) {
    object <- Seurat::DietSeurat(object)

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
    non_ensg_genes <- !grepl("^ENSG", features$gene_name)
    keep_genes <- keep_genes & non_mt_genes & non_ensg_genes

    # Update gene names and filter duplicates
    filtered_object <- object[keep_genes, ]

    assayobj <- filtered_object@assays$RNA

    # List of new gene names
    new_names <- features$gene_name[keep_genes]

    # Assign new names to the features slot of the temporary object
    rownames(assayobj@meta.features) <- new_names

    matrix_n <- SeuratObject::GetAssayData(assayobj)
    matrix_n@Dimnames[[1]] <- new_names

    # Assign new gene names to the temporary object
    assayobj@counts <- matrix_n
    assayobj@data <- matrix_n

    # Overwrite the RNA assay of the original object with the updated assay
    object@assays[["RNA"]] <- assayobj

    object <- object[rowSums(object) > 0, ]

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


#### 5. Plot gene expression scatter plot for each tissue ####

# Connect to Ensembl
ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
# Function to retrieve gene biotypes and prepare data for plotting
get_gene_biotypes <- function(object_filter) {
    gene_totals <- rowSums(object_filter)
    df <- data.frame(
        gene = rownames(object_filter),
        sum_exp = gene_totals,
        stringsAsFactors = FALSE
    )
    df$ENSG <- grepl("^ENSG", df$gene)

    # Connect to Ensembl
    ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

    # Get gene biotypes for non-ensg genes
    gene_annotations <- getBM(
        attributes = c("external_gene_name", "gene_biotype"),
        filters = "external_gene_name",
        values = df$gene[!df$ENSG],
        mart = ensembl
    )

    # Merge annotation with our data frame
    df <- left_join(df, gene_annotations, by = c("gene" = "external_gene_name"))

    # For ENSG genes, assign a custom label (here, "ENSG")
    df$gene_biotype <- ifelse(df$ENSG, "ENSG", df$gene_biotype)
    df$gene_biotype[is.na(df$gene_biotype)] <- "other" # in case any remain NA

    # Create a grouping variable for ENSG genes (custom label), protein_coding, lncRNA, other (any gene not falling in above)
    df$major_group <- ifelse(df$gene_biotype %in% c("protein_coding", "lncRNA"),
        df$gene_biotype, df$gene_biotype
    )
    # For clarity, ENSG remains "ENSG" and any remaining labels become "other"
    df$major_group[!(df$major_group %in% c("ENSG", "protein_coding", "lncRNA"))] <- "other"

    return(df)
}

# Function to generate gene expression scatter plot for a given tissue
generate_gene_expression_plot <- function(df, tissue_name) {
    # Define colors for each group (choose colors as you prefer)
    group_colors <- c("ENSG" = "red", "protein_coding" = "blue", "lncRNA" = "green", "other" = "black")

    # Plot the gene totals, coloring by the group
    p <- ggplot(df, aes(x = gene, y = sum_exp, color = major_group)) +
        geom_point() +
        scale_color_manual(values = group_colors) +
        labs(
            title = paste("Total Gene Expression Across All Cells -", tissue_name),
            x = "Gene",
            y = "Sum of expression"
        ) +
        theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())

    return(p)
}

# Retrieve gene biotypes and save to a named list
gene_df <- list(
    blood = get_gene_biotypes(blood_object_filter),
    lung = get_gene_biotypes(lung_object_filter),
    fat = get_gene_biotypes(fat_object_filter),
    kidney = get_gene_biotypes(kidney_object_filter),
    liver = get_gene_biotypes(liver_object_filter)
)

# Generate plots for all tissues
gex_plots <- list(
    blood = generate_gene_expression_plot(gene_df$blood, "blood"),
    lung = generate_gene_expression_plot(gene_df$lung, "lung"),
    fat = generate_gene_expression_plot(gene_df$fat, "fat"),
    kidney = generate_gene_expression_plot(gene_df$kidney, "kidney"),
    liver = generate_gene_expression_plot(gene_df$liver, "liver")
)

# Save all expression plots to a single PDF
pdf("plots/all_tissues_raw_exp_scatter.pdf", width = 10, height = 5)
for (plot in gex_plots) {
    print(plot)
}
dev.off()


#### 6. Plot batch effect ####
inspect_batch_effect <- function(object, tissue) {
    object <- Seurat::FindVariableFeatures(object, verbose = FALSE)
    object <- Seurat::ScaleData(object, verbose = FALSE)
    object <- Seurat::RunPCA(object, npcs = 30, verbose = FALSE, reduction.name = "pca")
    object <- Seurat::RunUMAP(object, reduction = "pca", dims = 1:20, reduction.name = "umap")
    object <- harmony::RunHarmony(object, group.by.vars = "donor_id", reduction = "umap")

    # Plot UMAP with batch effect
    p1 <- Seurat::DimPlot(object, reduction = "umap", group.by = "donor_id", pt.size = 0.75) +
        theme_minimal() +
        ggtitle(paste("Before Batch Correction - "), tissue)
    p2 <- Seurat::DimPlot(object, reduction = "harmony", group.by = "donor_id", pt.size = 0.75) +
        theme_minimal() +
        ggtitle(paste("After Batch Correction - "), tissue)
    donor_plot <- p1 + p2

    p3 <- Seurat::DimPlot(object, reduction = "umap", group.by = "cell_type", pt.size = 0.75) +
        theme_minimal() +
        ggtitle(paste("Before Batch Correction - "), tissue)
    p4 <- Seurat::DimPlot(object, reduction = "harmony", group.by = "cell_type", pt.size = 0.75) +
        theme_minimal() +
        ggtitle(paste("After Batch Correction - "), tissue)
    celltype_plot <- p3 + p4

    return(list(donor_plot = donor_plot, celltype_plot = celltype_plot))
}

batch_plots <- list(
    blood = inspect_batch_effect(blood_object_filter, "blood"),
    lung = inspect_batch_effect(lung_object_filter, "lung"),
    fat = inspect_batch_effect(fat_object_filter, "fat"),
    # kidney = inspect_batch_effect(kidney_object_filter, "kidney"),
    liver = inspect_batch_effect(liver_object_filter, "liver")
)

# Save all batch effect plots to a single PDF
pdf("plots/all_tissues_batch_effect_umap.pdf", width = 10, height = 5)
for (plot in batch_plots) {
    print(plot)
}
dev.off()




# amount of genes with non-zero expression across cells (tissue-specific)
# > length(rownames(blood_object_filter[rowSums(blood_object_filter) > 0, ]))
# [1] 56055
# > length(rownames(lung_object_filter[rowSums(lung_object_filter) > 0, ]))
# [1] 51646
# > length(rownames(liver_object_filter[rowSums(liver_object_filter) > 0, ]))
# [1] 43511
# > length(rownames(fat_object_filter[rowSums(fat_object_filter) > 0, ]))
# [1] 45310
# > length(rownames(kidney_object_filter[rowSums(kidney_object_filter) > 0, ]))
# [1] 45654
