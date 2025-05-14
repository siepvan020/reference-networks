#!/usr/bin/env Rscript


# Install required packages
if (!requireNamespace("ggplot2", quietly = TRUE)) install.packages("ggplot2")
if (!requireNamespace("Matrix", quietly = TRUE)) install.packages("Matrix")

if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
if (!requireNamespace("msigdbr", quietly = TRUE)) BiocManager::install("msigdbr")
if (!requireNamespace("fgsea", quietly = TRUE)) BiocManager::install("fgsea")

# Load required packages
library(ggplot2)
library(ggrepel)
library(patchwork)
library(msigdbr)
library(fgsea)
library(Matrix)
library(dplyr)
library(stringr)
library(Seurat)
library(uwot)
library(tidyr)
library(biomaRt)

# Set working directory
setwd("/div/pythagoras/u1/siepv/siep/Pipeline")

#### 1. Load and format data ####
load("output/networks/final/bloodScorpionOutput.Rdata")
load("output/networks/final/lungScorpionOutput.Rdata")
load("output/networks/final/fatScorpionOutput.Rdata")
load("output/networks/final/kidneyScorpionOutput.Rdata")
load("output/networks/final/liverScorpionOutput.Rdata")
load("output/networks/final/s_intestineScorpionOutput.Rdata")
load("output/networks/final/l_intestineScorpionOutput.Rdata")


#### 2. Calculate indegrees & outdegrees ####

calculate_indegrees <- function(output, tissue_name) {
    indegrees <- data.frame(gene = colnames(output[[1]]))
    for (celltype in names(output)) {
        regnet <- output[[celltype]]
        indegrees[[paste0(tissue_name, "_", celltype)]] <- colSums(regnet)
    }
    return(indegrees)
}

blood_indegrees <- calculate_indegrees(blood_output, "blood")
lung_indegrees <- calculate_indegrees(lung_output, "lung")
fat_indegrees <- calculate_indegrees(fat_output, "fat")
kidney_indegrees <- calculate_indegrees(kidney_output, "kidney")
liver_indegrees <- calculate_indegrees(liver_output, "liver")
s_intestine_indegrees <- calculate_indegrees(s_intestine_output, "s_intestine")
l_intestine_indegrees <- calculate_indegrees(l_intestine_output, "l_intestine")

# Combine blood and lung indegrees and outdegrees into one dataframe
combined_indegrees <- Reduce(function(x, y) merge(x, y, by = "gene", all = TRUE), list(blood_indegrees, lung_indegrees, fat_indegrees, kidney_indegrees, liver_indegrees, s_intestine_indegrees, l_intestine_indegrees))
combined_indegrees <- combined_indegrees[, c("gene", sort(colnames(combined_indegrees)[-1]))]

# Generate PCA and UMAP and plot only UMAP on different conditions
perform_dimensionality_reduction <- function(data, output_prefix, analysis_type, avgdepth = NULL, n_neighbors = 8, min_dist = 0.1) {
    # Prepare matrix
    mat                 <- as.matrix(data[, -1])
    rownames(mat)       <- data$gene
    mat[is.na(mat)]     <- 0
    mat                 <- t(mat)
    mat                 <- scale(mat)

    # Perform PCA
    pca_res             <- prcomp(mat, center = TRUE, scale. = TRUE)
    pca_df              <- as.data.frame(pca_res$x)
    pca_df$condition    <- rownames(pca_df)
    pca_df$tissue       <- ifelse(grepl("^s_intestine", pca_df$condition), "s_intestine", 
                                ifelse(grepl("^l_intestine", pca_df$condition), "l_intestine", 
                                sub("_.*$", "", pca_df$condition)))
    pca_df$celltype     <- sub("^s_intestine_|^l_intestine_|^[^_]+_", "", pca_df$condition)
    if (!is.null(avgdepth)) { pca_df$counts <- avgdepth$celltype_exp }
    umap_input          <- pca_df[, grep("^PC", colnames(pca_df))]

    # Perform UMAP
    set.seed(42)
    umap_res            <- uwot::umap(umap_input, n_neighbors = n_neighbors, min_dist = min_dist)
    umap_df             <- as.data.frame(umap_res)
    colnames(umap_df)   <- c("UMAP1", "UMAP2")
    umap_df$condition   <- rownames(umap_df)
    umap_df$tissue      <- pca_df$tissue
    umap_df$celltype    <- pca_df$celltype
    if (!is.null(avgdepth)) { umap_df$counts <- avgdepth$celltype_exp }

    save_umap <- function(color, shape = NULL, use_gradient = FALSE) {
        mapping <- aes(UMAP1, UMAP2, color = !!rlang::sym(color))
        if (!use_gradient && !is.null(shape)) {
            mapping <- modifyList(mapping, aes(shape = !!rlang::sym(shape)))
        }

        plt <- ggplot(umap_df, mapping) +
            geom_point(size = 4) +
            labs(
                title = paste("UMAP of", analysis_type),
                x = "UMAP1", y = "UMAP2",
                colour = color,
                shape = if (!use_gradient && !is.null(shape)) shape else NULL
            ) +
            theme_minimal()

        if (use_gradient) {
            plt <- plt +
                scale_colour_gradient(low = "white", high = "blue") +
                guides(color = guide_colorbar(title = "Mean UMI per cell type"))
            suffix <- "_counts"
        } else {
            plt <- plt + scale_colour_discrete() +
                scale_shape_manual(values = c(16, 17, 15, 18, 7, 8, 4))
            suffix <- ""
        }

        ggsave(paste0(
            "output/analysis/dimreduction/umap/",
            output_prefix, "_umap", suffix, ".pdf"
        ),
        plt,
        width = 12, height = 6
        )
    }

    save_umap("celltype", "tissue")
    if (!is.null(avgdepth)) {
        save_umap("counts", use_gradient = TRUE)
    }

    return(list(pca = pca_df, umap = umap_df, mat = mat, pca_res = pca_res))
}

indegree_reducs <- perform_dimensionality_reduction(combined_indegrees, "indegree", "Indegrees") # Only plot UMAP colored and shaped by celltype and tissue

indegree_reducs_cnts <- perform_dimensionality_reduction(combined_indegrees, "indegree", "Indegrees", avgdepth) # run this after the expression data is calculated
exp_reducs <- perform_dimensionality_reduction(combined_exp_matrix, "expression", "Expression", avgdepth) # run this after the expression data is calculated

# Calculate correlation matrix
cor_indegree_matrix <- cor(combined_indegrees[, -1], use = "pairwise.complete.obs")

# Convert the correlation matrix to a long format for ggplot
cor_indegree_long <- as.data.frame(as.table(cor_indegree_matrix))

# Function to plot the correlation matrix with values in the tiles
plot_correlation_matrix <- function(data, title = "Correlation Matrix", gradient_limits = c(NA, NA), labels = FALSE) {
    plt <- ggplot(data, aes(Var1, Var2, fill = Freq)) +
        geom_tile() +
        scale_fill_gradient(
            low = "white",
            high = "steelblue",
            name = "Pearson's r value",
            limits = if (all(is.na(gradient_limits))) NULL else gradient_limits
        ) +
        theme_minimal() +
        theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
        labs(x = "Cell Types", y = "Cell Types", title = title)
    if (labels) {
        plt <- plt + geom_text(aes(label = round(Freq, 3)), color = "black", size = 3)
    }
    plt
}

cor_indegree_plot <- plot_correlation_matrix(cor_indegree_long,
    title = "Correlation Matrix of Indegrees"
)
ggsave("output/analysis/correlation/correlation_matrix_indegrees.pdf", cor_indegree_plot, width = 12, height = 6)


#### 3. Calculate linear regression and residuals between conditions ####

# Define comparisons
comparisons <- list(
    list(name = "cd4_blood_vs_s_intestine", x = "blood_cd4_positive_t_cell", y = "s_intestine_cd4_positive_t_cell"),
    list(name = "monocyte_blood_vs_liver", x = "blood_monocyte", y = "liver_monocyte")
)

# Function to fit linear model and extract + rank residuals descending
calculate_residuals <- function(data, x_con, y_con) {
    fit <- lm(data[[y_con]] ~ data[[x_con]], data = data) # Fit linear model
    res <- residuals(fit)
    names(res) <- data$gene
    return(list(fit = fit, residuals = res[order(-res)])) # Return fit and ranked residuals
}

plot_linear_regression <- function(data, x_con, y_con, comp_name, residuals, fit, analysis) {
    # Compute prediction intervals
    data$fit <- predict(fit, newdata = data)

    # Identify top 5 and bottom 5 genes
    extreme_genes <- c(head(names(residuals), 5), tail(names(residuals), 5))

    # Assign colors and labels for extreme genes
    data$color <- ifelse(data$gene %in% extreme_genes, ifelse(data$gene %in% head(names(residuals), 5), "red", "blue"), "#2c2c2cc4")
    data$label <- ifelse(data$gene %in% extreme_genes, data$gene, NA)

    # Build the plot
    p <- ggplot(data, aes(x = x_con, y = y_con)) +
        geom_point(alpha = 0.5, aes(color = color)) +
        geom_point(data = filter(data, gene %in% extreme_genes), aes(color = color), alpha = 0.75, size = 1.25) +
        geom_line(aes(y = fit), color = "red", linewidth = 0.75) +
        geom_density_2d(color = "gray", alpha = 0.5) +
        geom_text_repel(
            data = filter(data, gene %in% extreme_genes),
            aes(label = label, color = color),
            fontface = 3,
            min.segment.length = 0,
            size = 3,
            segment.size = 0.1,
            max.overlaps = 75,
            bg.color = "white"
        ) +
        scale_color_identity() +
        labs(
            title = paste("Linear Regression", analysis, ":", comp_name),
            x = paste(x_con, analysis),
            y = paste(y_con, analysis)
        ) +
        theme_light() +
        theme(legend.position = "none")
}

# Calculate residuals for each comparison
residuals_list <- list()
fit_list <- list()
plot_list <- list()
for (comp in comparisons) {
    cat("Calculating residuals for:", comp$name, "\n")
    data <- combined_indegrees[complete.cases(combined_indegrees[[comp$x]], combined_indegrees[[comp$y]]), ] # Only keep complete cases between conditions
    result <- calculate_residuals(data, comp$x, comp$y)
    residuals_list[[comp$name]] <- result$residuals
    fit_list[[comp$name]] <- result$fit
    plot_list[[comp$name]] <- plot_linear_regression(data, comp$x, comp$y, comp$name, result$residuals, result$fit, "indegree")
}

# Save plots of the linear regression models
ggsave("output/analysis/linear_regression/indegree/cd4_blood_vs_s_intestine.pdf", plot_list[[2]], width = 12, height = 6)
ggsave("output/analysis/linear_regression/indegree/monocyte_blood_vs_liver.pdf", plot_list[[4]], width = 12, height = 6)



#### 4. Run GSEA on ranked residuals ####

run_gsea <- function(residuals_list, minSize = 10, maxSize = 500, nproc = 1) {
    gsea_results <- list()
    options(warn = 1)

    # Load gene sets for GO:BP and GO:MF using gene symbols
    msigdb_BP <- msigdbr(species = "Homo sapiens", collection = "C5", subcollection = "BP")[, c("gs_name", "gene_symbol")]
    msigdb_MF <- msigdbr(species = "Homo sapiens", collection = "C5", subcollection = "MF")[, c("gs_name", "gene_symbol")]

    # Convert to fgsea-compatible format (list of gene sets)
    msigdb_BP <- split(msigdb_BP$gene_symbol, msigdb_BP$gs_name)
    msigdb_MF <- split(msigdb_MF$gene_symbol, msigdb_MF$gs_name)

    set.seed(42)
    for (comp in names(residuals_list)) {
        cat("Running GSEA for:", comp, "\n")

        ranked_genes <- residuals_list[[comp]] # Get ranked gene list

        # Run fgsea for GO:BP
        gsea_results[[paste0("BP_", comp)]] <- fgsea(
            pathways = msigdb_BP,
            stats = ranked_genes,
            minSize = minSize,
            maxSize = maxSize,
            nproc = nproc
        )

        # Run fgsea for GO:MF
        gsea_results[[paste0("MF_", comp)]] <- fgsea(
            pathways = msigdb_MF,
            stats = ranked_genes,
            minSize = minSize,
            maxSize = maxSize,
            nproc = nproc
        )
    }
    return(gsea_results)
}

gsea_results_indegree <- run_gsea(residuals_list)

# Function to get top up and down pathways
get_top_pathways <- function(gsea_result) {
    sig_pathways <- gsea_result %>%
        filter(padj < 0.05)

    topPathwaysUp <- sig_pathways %>%
        filter(NES > 0) %>%
        arrange(desc(NES)) %>%
        head(10) %>%
        mutate(pathway = factor(pathway, levels = pathway))

    topPathwaysDown <- sig_pathways %>%
        filter(NES < 0) %>%
        arrange(NES) %>%
        head(10) %>%
        mutate(pathway = factor(pathway, levels = pathway))

    return(list(up = topPathwaysUp, down = topPathwaysDown, all = sig_pathways))
}

# Function to plot top pathways
plot_top_pathways <- function(topPathways, comparison_name) {
    for (dir in c("up", "down")) {
        topPathways[[dir]] <- topPathways[[dir]] %>%
            mutate(
                pathway = gsub("_", " ", as.character(pathway)),
                pathway = str_wrap(pathway, width = 50, whitespace_only = FALSE), # Wrap long pathway names
                pathway = gsub(" ", "_", pathway)
            )
    }

    p1 <- ggplot(topPathways$up, aes(x = NES, y = reorder(pathway, NES), size = size, fill = NES)) +
        geom_point(shape = 21) +
        scale_size_area(max_size = 10) +
        scale_fill_gradient(low = "white", high = "red") +
        theme(
            axis.text.y = element_text(size = 16),
            axis.title.y = element_blank(),
            plot.title = element_text(hjust = 1, face = "bold")
        ) +
        xlim(min(topPathways$up$NES), max(topPathways$up$NES)) +
        xlab("NES") +
        ggtitle(paste("Top 10 positively enriched pathways in", comparison_name))

    p2 <- ggplot(topPathways$down, aes(x = NES, y = reorder(pathway, NES), size = size, fill = NES)) +
        geom_point(shape = 21) +
        scale_size_area(max_size = 10) +
        scale_fill_gradient(low = "blue", high = "white") +
        theme(
            axis.text.y = element_text(size = 16),
            axis.title.y = element_blank(),
            plot.title = element_text(hjust = 1, face = "bold")
        ) +
        xlim(min(topPathways$down$NES), max(topPathways$down$NES)) +
        xlab("NES") +
        ggtitle(paste("Top 10 negatively enriched pathways in", comparison_name))

    return(list(up = p1, down = p2))
}

# Generate plots for all comparisons in gsea_results and save to a list
gsea_plots <- list()
top_pathways_list <- list()
iteration <- 1
for (comp in names(gsea_results_indegree)) {
    # Generate plots
    cat("Plotting GSEA results for:", comp, "\n")
    topPathways <- get_top_pathways(gsea_results_indegree[[comp]]) # Returns up, down and all significant pathways
    gsea_plots[[comp]] <- plot_top_pathways(topPathways, comp)
    top_pathways_list[[comp]] <- topPathways

    # Save plots
    updotplot <- gsea_plots[[comp]]$up
    downdotplot <- gsea_plots[[comp]]$down
    pdf(file = paste0("output/analysis/gsea_dotplots/indegree/", iteration, "_", comp, "_dotplot.pdf"), width = 20, height = 13)
    print(patchwork::wrap_plots(updotplot, downdotplot, ncol = 1))
    dev.off()
    iteration <- iteration + 1
}


#### 5. Run analysis on expression data ####

blood_object <- readRDS("output/preprocessing/blood_prepped.rds")
lung_object <- readRDS("output/preprocessing/lung_prepped.rds")
fat_object <- readRDS("output/preprocessing/fat_prepped.rds")
kidney_object <- readRDS("output/preprocessing/kidney_prepped.rds")
liver_object <- readRDS("output/preprocessing/liver_prepped.rds")
s_intestine_object <- readRDS("output/preprocessing/s_intestine_prepped.rds")
l_intestine_object <- readRDS("output/preprocessing/l_intestine_prepped.rds")


calculate_celltype_depth <- function(object, tissue) {
    cond_depths <- sapply(levels(Idents(object)), function(ct) {
        cells <- WhichCells(object, idents = ct)
        mean(Matrix::colSums(object[["RNA"]]@counts[, cells]))
    })
    names(cond_depths) <- paste0(tissue, "_", names(cond_depths))
    return(cond_depths)
}

# Per tissue
blood_depths <- calculate_celltype_depth(blood_object, "blood")
lung_depths <- calculate_celltype_depth(lung_object, "lung")
fat_depths <- calculate_celltype_depth(fat_object, "fat")
kidney_depths <- calculate_celltype_depth(kidney_object, "kidney")
liver_depths <- calculate_celltype_depth(liver_object, "liver")
s_intestine_depths <- calculate_celltype_depth(s_intestine_object, "s_intestine")
l_intestine_depths <- calculate_celltype_depth(l_intestine_object, "l_intestine")
# Combine into one vector
all_depths <- c(
    blood_depths, lung_depths, fat_depths,
    kidney_depths, liver_depths,
    s_intestine_depths, l_intestine_depths
)
avgdepth <- data.frame(celltype_exp = all_depths)

#
#
#
##### Calculate sum of expression for each gene in each cell type - pseudobulks #####
calculate_expression_sum <- function(object, tissue) {
    sum_df <- data.frame(gene = rownames(object[["RNA"]]@counts))
    for (ct in levels(Idents(object))) {
        ct_subset <- subset(object, idents = ct)
        sum <- rowSums(ct_subset[["RNA"]]@counts)
        sum_df[[paste0(tissue, "_", ct)]] <- sum
    }
    sum_df
}

blood_sum_exp <- calculate_expression_sum(blood_object, "blood")
lung_sum_exp <- calculate_expression_sum(lung_object, "lung")
fat_sum_exp <- calculate_expression_sum(fat_object, "fat")
kidney_sum_exp <- calculate_expression_sum(kidney_object, "kidney")
liver_sum_exp <- calculate_expression_sum(liver_object, "liver")
s_intestine_sum_exp <- calculate_expression_sum(s_intestine_object, "s_intestine")
l_intestine_sum_exp <- calculate_expression_sum(l_intestine_object, "l_intestine")

# Combine all data frames into one big data frame
combined_exp_matrix <- Reduce(
    function(x, y) merge(x, y, by = "gene", all = TRUE),
    list(
        blood_sum_exp, lung_sum_exp, fat_sum_exp,
        kidney_sum_exp, liver_sum_exp,
        s_intestine_sum_exp, l_intestine_sum_exp
    )
)
combined_exp_matrix <- combined_exp_matrix[, c("gene", sort(colnames(combined_exp_matrix)[-1]))]

cor_exp_matrix <- cor(combined_exp_matrix[, -1], use = "pairwise.complete.obs")


cor_exp_long <- as.data.frame(as.table(cor_exp_matrix))
difference_long <- as.data.frame(as.table(cor_indegree_matrix - cor_exp_matrix))

cor_exp_plot <- plot_correlation_matrix(cor_exp_long, title = "Correlation Matrix of Expression", gradient_limits = c(0, 1))
ggsave("output/analysis/correlation/correlation_matrix_expression.pdf", cor_exp_plot, width = 12, height = 6)
cor_dif_plot <- plot_correlation_matrix(difference_long, title = "Difference in Correlation Matrix of Indegrees and Expression", gradient_limits = c(min(difference_long$Freq), 1))
ggsave("output/analysis/correlation/correlation_matrices_difference.pdf", cor_dif_plot, width = 12, height = 6)


##### Calculate residuals for each comparison #####
gex_residuals <- list()
gex_fit <- list()
gex_plot <- list()
for (comp in comparisons) {
    cat("Calculating gex residuals for:", comp$name, "\n")
    data <- combined_exp_matrix[complete.cases(combined_exp_matrix[[comp$x]], combined_exp_matrix[[comp$y]]), ] # Only keep complete cases between conditions
    result <- calculate_residuals(data, comp$x, comp$y)
    gex_residuals[[comp$name]] <- result$residuals
    gex_fit[[comp$name]] <- result$fit
    gex_plot[[comp$name]] <- plot_linear_regression(data, comp$x, comp$y, comp$name, result$residuals, result$fit, "sum of expression per gene")
}

# Save plots of the linear regression models
ggsave("output/analysis/linear_regression/gex/cd4_blood_vs_s_intestine.pdf", gex_plot[[1]], width = 12, height = 6)
ggsave("output/analysis/linear_regression/gex/monocyte_blood_vs_liver.pdf", gex_plot[[2]], width = 12, height = 6)

# Run GSEA on the expression data
gsea_results_gex <- run_gsea(gex_residuals)

# Generate plots for all comparisons in gsea_results and save to a list
gsea_gex_plots <- list()
gex_top_pathways_list <- list()
iteration <- 1
for (comp in names(gsea_results_gex)) {
    # Generate plots
    cat("Plotting GSEA results for:", comp, "\n")
    topPathways <- get_top_pathways(gsea_results_gex[[comp]])
    gsea_gex_plots[[comp]] <- plot_top_pathways(topPathways, comp)
    gex_top_pathways_list[[comp]] <- topPathways

    # Save plots
    updotplot <- gsea_gex_plots[[comp]]$up
    downdotplot <- gsea_gex_plots[[comp]]$down
    pdf(file = paste0("output/analysis/gsea_dotplots/gex/", iteration, "_", comp, "_dotplot.pdf"), width = 17, height = 15)
    print(patchwork::wrap_plots(updotplot, downdotplot, ncol = 1))
    dev.off()
    iteration <- iteration + 1
}
