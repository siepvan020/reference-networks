#!/usr/bin/env Rscript


#### 1. Setup ####

progress_file <- "/div/pythagoras/u1/siepv/siep/Analysis_v2/output/log/preprocess.log"

# Start fresh log or append
cat("- Starting script at -", format(Sys.time()), "\n", file = progress_file)

if (!requireNamespace("SeuratObject", quietly = TRUE)) remotes::install_version("SeuratObject", "4.1.4", repos = c("https://satijalab.r-universe.dev", getOption("repos")))
if (!requireNamespace("Seurat", quietly = TRUE)) remotes::install_version("Seurat", "4.4.0", repos = c("https://satijalab.r-universe.dev", getOption("repos")))

library(Seurat)
library(SeuratObject)
library(biomaRt)
library(dplyr)
library(ggplot2)
library(sva)

library(yaml)
library(glue)
library(purrr)

# Set working directory
setwd("/div/pythagoras/u1/siepv/siep/Analysis_v2")

# Load config file
config <- yaml::read_yaml("data/config/config_full.yaml")


#### 2. Load data & filter cells ####
# Convert NULL to empty character vector
.null_to_chr <- function(x) if (is.null(x)) character(0) else x

# Reads one tissue object, merges and corrects celltype annotations, then
# keeps only the desired cell types and returns the filtered Seurat object.
process_tissue <- function(spec, tissue) {
    tryCatch({
        object <- readRDS(glue("data/rds/{spec$file}"))

        mappings <- map(spec$cells, .null_to_chr) |> # ensure vectors
            keep(~ length(.x) > 0) # drop empty ones

        object <- merge_cell_types(object, mappings) # Returns object with updated cell types in metadata
        object <- object[, object@meta.data$cell_type %in% names(spec$cells)]
        object
    }, error = function(e) {
        cat(glue("Error processing tissue '{tissue}': {e$message}"), format(Sys.time()), "\n", file = progress_file)
        NULL
    })
}

# Function to merge cell types into a new label based on the config mapping
merge_cell_types <- function(object, cell_type_mappings) {
    object <- Seurat::DietSeurat(object)

    object@meta.data$cell_type <- as.character(object@meta.data$cell_type)

    for (new_label in names(cell_type_mappings)) {
        old_labels <- cell_type_mappings[[new_label]]
        object@meta.data$cell_type[object@meta.data$cell_type %in% old_labels] <- new_label
    }
    object@meta.data$cell_type <- as.factor(object@meta.data$cell_type)
    Idents(object) <- "cell_type"

    object
}

# Call the process_tissue function for each tissue in the config
objects <- purrr::imap(config, process_tissue)
cat("- Loaded and cell-filtered Seurat objects -", format(Sys.time()), "\n", file = progress_file, append = TRUE)


#### 3. Filter genes ####
# Function to update gene names and filter the expression matrix
update_gene_names_and_filter <- function(object, keep_genes, features) {
    non_mt <- !object@assays$RNA@meta.features$mt # To get correct sum of expression PLOT, comment out this line
    non_ensg <- !grepl("^ENSG", features$gene_name) # To get correct sum of expression PLOT, comment out this line
    keep_vec <- keep_genes & non_mt & non_ensg # To get correct sum of expression PLOT, comment out this line

    obj_filt <- object[keep_vec, ]
    assay <- obj_filt@assays$RNA
    new_names <- features$gene_name[keep_vec]

    rownames(assay@meta.features) <- new_names
    mat <- SeuratObject::GetAssayData(assay)
    mat@Dimnames[[1]] <- new_names
    assay@counts <- mat
    assay@data <- mat
    obj_filt@assays[["RNA"]] <- assay

    obj_filt <- obj_filt[rowSums(obj_filt) > 0, ]
    obj_filt
}

# Extract feature names and Ensembl IDs from the first object in the list,
# assuming every object has the same features
template_obj <- objects[[1]]
features <- data.frame(
    gene_name  = sub("_ENSG[0-9]+$", "", template_obj@assays$RNA@meta.features$feature_name),
    ensembl_id = rownames(template_obj)
)

# Identify duplicate gene names as not all ENSG IDs are unique
dup_ids <- features$gene_name[duplicated(features$gene_name) | duplicated(features$gene_name, fromLast = TRUE)]
genes_to_exclude <- unique(features$ensembl_id[features$gene_name %in% dup_ids])
keep_genes <- !rownames(template_obj) %in% genes_to_exclude

# Apply the gene filtering function to each object in the objects list
objects <- purrr::map(objects, update_gene_names_and_filter, keep_genes, features)
cat("- Gene-filtered Seurat objects -", format(Sys.time()), "\n", file = progress_file, append = TRUE)

# Function to generate UMAP plot to visualize the batch effect
inspect_batch_effect <- function(object, tissue, stage = c("before", "after")) {
    stage <- match.arg(stage)

    # Set reduction name based on stage
    reduction_name <- if (stage == "before") "umap_before" else "umap_after"

    # Run Seurat steps
    object <- Seurat::FindVariableFeatures(object, verbose = FALSE)
    object <- Seurat::ScaleData(object, verbose = FALSE)
    object <- Seurat::RunPCA(object, npcs = 30, verbose = FALSE, reduction.name = paste0("pca_", stage))
    object <- Seurat::RunUMAP(object, reduction = paste0("pca_", stage), dims = 1:20, reduction.name = reduction_name)

    # Plot UMAP grouped by donor_id
    plot_donor <- Seurat::DimPlot(object, reduction = reduction_name, group.by = "donor_id", pt.size = 0.75) +
        theme_minimal() +
        ggtitle(paste(toupper(stage), "Batch Correction -", tissue, "- Donor ID"))

    # Plot UMAP grouped by cell_type
    plot_cell_type <- Seurat::DimPlot(object, reduction = reduction_name, group.by = "cell_type", pt.size = 0.75) +
        theme_minimal() +
        ggtitle(paste(toupper(stage), "Batch Correction -", tissue, "- Cell Type"))

    return(list(plot_donor = plot_donor, plot_cell_type = plot_cell_type))
}


#### 4. Create UMAP plots BEFORE batch correction ####
batch_plots_before <- purrr::imap(
    objects,
    function(.x, .y) {
        if (length(unique(.x$donor_id)) > 1) {
            inspect_batch_effect(.x, .y, stage = "before")
        } else {
            NULL
        }
    }
)


#### 5. Batch correction with ComBat_seq and save objects ####
perform_save_batch_correction <- function(obj, tissue) {
    tryCatch({
        donors <- obj$donor_id
        if (length(unique(donors)) > 1) { # Skip tissues with one donor
            cat(glue("   - Batch correcting {tissue} object"), format(Sys.time()), "\n", file = progress_file, append = TRUE)
            counts_corr <- sva::ComBat_seq(as.matrix(obj@assays$RNA@counts),
                batch = donors
            )
            obj[["RNA"]]@counts <- counts_corr
        }
        saveRDS(obj, glue::glue("output/preprocessing/all/{tissue}_prepped.rds"))
        cat(glue("   - Saved batch-corrected object for {tissue}"), format(Sys.time()), "\n", file = progress_file, append = TRUE)
        obj
    }, error = function(e) {
        cat(glue("Error processing batch correction for tissue '{tissue}': {e$message}"), format(Sys.time()), "\n", file = progress_file, append = TRUE)
        NULL
    })
}

cat("- Starting batch correction and saving using ComBat_seq -", format(Sys.time()), "\n", file = progress_file, append = TRUE)
objects <- purrr::imap(
    objects,
    ~ perform_save_batch_correction(.x, .y)
)
cat("- Completed batch correction and saving -", format(Sys.time()), "\n", file = progress_file, append = TRUE)


#### 7. Create UMAP plots AFTER batch correction ####
batch_plots_after <- purrr::imap(
    objects,
    function(.x, .y) {
        if (length(unique(.x$donor_id)) > 1) {
            inspect_batch_effect(.x, .y, stage = "after")
        } else {
            NULL
        }
    }
)


#### 8. Save UMAP plots to PDF ####

# Combine before and after batch effect plots for each tissue without a loop
combined_batch_plots <- purrr::imap(
    batch_plots_before,
    ~ list(
        combined_donor_plot = .x$plot_donor + batch_plots_after[[.y]]$plot_donor,
        combined_cell_type_plot = .x$plot_cell_type + batch_plots_after[[.y]]$plot_cell_type
    )
)

# Save combined plots to PDF
pdf("output/preprocessing/plots/all_tissues_batch_effect_umap.pdf", width = 10, height = 5)
purrr::iwalk(combined_batch_plots, function(plots) {
    print(plots$combined_donor_plot)
    print(plots$combined_cell_type_plot)
})
dev.off()


#### 9. Generate gene expression scatter plot ####

# Connect to Ensembl
# ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
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
    # Define colors for each group
    group_colors <- c("ENSG" = "red", "protein_coding" = "blue", "lncRNA" = "green", "other" = "black")

    # Plot the gene totals, coloring by the group
    p <- ggplot(df, aes(x = gene, y = sum_exp, color = major_group)) +
        geom_point() +
        scale_color_manual(values = group_colors) +
        labs(
            title = paste("Total gene expression across subset of cells in", tissue_name),
            x = "Gene A-Z",
            y = "Sum of expression",
            color = "Gene biotype"
        ) +
        theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())

    return(p)
}

# # Get gene biotypes and generate expression plots for each tissue
# gene_df <- purrr::imap(objects, ~ get_gene_biotypes(.x))
# gex_plots <- purrr::imap(gene_df, ~ generate_gene_expression_plot(.x, .y))

# # Save the gene expression scatter plots to a PDF
# pdf("output/preprocessing/plots/all_tissues_raw_exp_scatter.pdf",
#     width = 10,
#     height = 3
# )
# purrr::iwalk(gex_plots, ~ print(.x))
# dev.off()
