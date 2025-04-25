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
# library(biomaRt)

# Set working directory
setwd("/div/pythagoras/u1/siepv/siep/Analysis_v2")

#### 1. Load and format data ####
load("output/networks/final/bloodScorpionOutput3.Rdata")
load("output/networks/final/lungScorpionOutput3.Rdata")
# load("final/fatScorpionOutput3.Rdata")
# load("final/kidneyScorpionOutput3.Rdata")
# load("final/liverScorpionOutput3.Rdata")


#### 2. Calculate indegrees & outdegrees ####

calculate_indegrees <- function(output, tissue_name) {
    indegrees <- data.frame(gene = colnames(output[[1]]$regNet))
    for (celltype in names(output)) {
        regnet <- output[[celltype]]$regNet
        indegrees[[paste0(tissue_name, "_", celltype)]] <- colSums(regnet)
    }
    return(indegrees)
}

blood_indegrees <- calculate_indegrees(blood_output, "blood")
lung_indegrees <- calculate_indegrees(lung_output, "lung")
# fat_indegrees <- calculate_indegrees(fat_output, "fat")
# kidney_indegrees <- calculate_indegrees(kidney_output, "kidney")
# liver_indegrees <- calculate_indegrees(liver_output, "liver")

# Combine blood and lung indegrees and outdegrees into one dataframe
combined_indegrees <- Reduce(function(x, y) merge(x, y, by = "gene", all = TRUE), list(blood_indegrees, lung_indegrees)) # , fat_indegrees, kidney_indegrees, liver_indegrees))


# try-out function to make pca and umap plots of indegrees and expression
perform_dimensionality_reduction <- function(data, output_prefix, analysis_type, colsum=NULL, n_neighbors = 8, min_dist = 0.1) {
    # Prepare matrix
    mat <- as.matrix(data[complete.cases(data), -1])
    rownames(mat) <- data$gene[complete.cases(data)]
    mat <- t(mat)

    # Perform PCA
    pca_res <- prcomp(mat, center = TRUE, scale. = TRUE)
    pca_df <- as.data.frame(pca_res$x)
    pca_df$condition <- rownames(pca_df)
    pca_df$counts <- colsum$celltype_exp

    # Perform UMAP
    set.seed(42)
    umap_res <- umap(pca_df, n_neighbors = n_neighbors, min_dist = min_dist)
    umap_df <- as.data.frame(umap_res)
    colnames(umap_df) <- c("UMAP1", "UMAP2")
    umap_df$condition <- rownames(umap_df)
    umap_df$counts <- colsum$celltype_exp

    # Plot PCA - change pc1 and pc2 to get different plots
    pca_plot <- ggplot(pca_df, aes(x = PC1, y = PC2, label = condition)) +
        geom_point(shape = 21, size = 4)+#, aes(fill = log10(counts + 1))) +
        geom_text_repel(size = 4, max.overlaps = 10) +
        # scale_fill_gradient(low = "white", high = "blue", name = "Log10 Counts") +
        labs(x = "PC1", y = "PC2", title = paste("PCA of", analysis_type)) +
        theme_minimal()
    ggsave(paste0("output/analysis/dimreduction/pca/", output_prefix, "_pca12.pdf"), pca_plot, width = 12, height = 6)

    umap_plot <- ggplot(umap_df, aes(x = UMAP1, y = UMAP2, label = condition))+#, color = celltype)) +
        geom_point(shape = 21, size = 4, aes(fill = log10(counts + 1))) +
        geom_text_repel(size = 4, max.overlaps = 10) +
        scale_fill_gradient(low = "white", high = "blue", name = "Log10 Counts") +
        labs(x = "UMAP1", y = "UMAP2", title = paste("UMAP of", analysis_type)) +
        theme_minimal()
    ggsave(paste0("output/analysis/dimreduction/umap/", output_prefix, "_umap.pdf"), umap_plot, width = 12, height = 6)

    return(list(pca = pca_df, umap = umap_df, mat = mat, pca_res = pca_res))
}

colsum <- data.frame(celltype_exp = colSums(combined_exp_matrix[complete.cases(combined_exp_matrix), ][, -1])) # run this after the expression data is calculated

indegree_reducs <- perform_dimensionality_reduction(combined_indegrees, "indegree", "Indegrees", colsum) # run this after the expression data is calculated
exp_reducs <- perform_dimensionality_reduction(combined_exp_matrix, "expression", "Expression", colsum) # run this after the expression data is calculated

# Calculate correlation matrix
cor_indegree_matrix <- cor(combined_indegrees[, -1], use = "pairwise.complete.obs")

# Convert the correlation matrix to a long format for ggplot
cor_long <- as.data.frame(as.table(cor_indegree_matrix))

# Function to plot the correlation matrix with values in the tiles
plot_correlation_matrix <- function(data, title = "Correlation Matrix", gradient_limits = c(NA, NA)) {
    ggplot(data, aes(Var1, Var2, fill = Freq)) +
        geom_tile(na.rm = TRUE) +
        geom_text(aes(label = round(Freq, 3)), color = "black", size = 3, na.rm = TRUE) +
        scale_fill_gradient(
            low = "white",
            high = "steelblue",
            name = "Pearson's r value",
            limits = if (all(is.na(gradient_limits))) NULL else gradient_limits
        ) +
        theme_minimal() +
        theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
        labs(x = "Cell Types", y = "Cell Types", title = title)
}

cor_indegree_plot <- plot_correlation_matrix(cor_long, title = "Correlation Matrix of Indegrees")
ggsave("output/analysis/correlation/run3_correlation_matrix_indegrees", cor_indegree_plot, width = 12, height = 6)


#### 3. Calculate linear regression and residuals between conditions ####

# Define comparisons
comparisons <- list(
    list(name = "cd4_blood_vs_lung", x = "blood_cd4_positive_t_cell", y = "lung_cd4_positive_t_cell"),
    list(name = "cd8_blood_vs_lung", x = "blood_cd8_positive_t_cell", y = "lung_cd8_positive_t_cell"),
    list(name = "monocyte_blood_vs_lung", x = "blood_monocyte", y = "lung_monocyte"),
    list(name = "lung_cd4_vs_cd8", x = "lung_cd4_positive_t_cell", y = "lung_cd8_positive_t_cell"),
    list(name = "lung_cd4_vs_monocyte", x = "lung_cd4_positive_t_cell", y = "lung_monocyte"),
    list(name = "lung_cd8_vs_monocyte", x = "lung_cd8_positive_t_cell", y = "lung_monocyte"),
    list(name = "blood_platelet_vs_lung_2pneumo", x = "blood_platelet", y = "lung_type_ii_pneumocyte"),
    list(name = "blood_monocyte_vs_lung_2pneumo", x = "blood_monocyte", y = "lung_type_ii_pneumocyte")
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
    p <- ggplot(data, aes_string(x = x_con, y = y_con)) +
        geom_point(alpha = 0.5, aes(color = color)) +
        geom_point(data = filter(data, gene %in% extreme_genes), aes(color = color), alpha = 0.75, size = 1.25) +
        geom_line(aes(y = fit), color = "red", size = 0.75) +
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
ggsave("output/analysis/linear_regression/run3/immune_cells_between_tissues3.pdf", patchwork::wrap_plots(A = plot_list[[1]], B = plot_list[[2]], C = plot_list[[3]], design = "AABB\n#CC#"), width = 12, height = 6)
ggsave("output/analysis/linear_regression/run3/immune_cells_lung3.pdf", patchwork::wrap_plots(A = plot_list[[4]], B = plot_list[[5]], C = plot_list[[6]], design = "AABB\n#CC#"), width = 12, height = 6)
ggsave("output/analysis/linear_regression/run3/tissue_specific_linear3.pdf", plot_list[[7]], width = 12, height = 6)
ggsave("output/analysis/linear_regression/run3/blood_monocyte_vs_lung_2pneumo.pdf", plot_list[[8]], width = 12, height = 6)


#### 4. Run GSEA on ranked residuals ####

run_gsea <- function(residuals_list, minSize = 10, maxSize = 500, nproc = 1) {
    gsea_results <- list()
    options(warn = 1)

    # Load gene sets for GO:BP and GO:MF using gene symbols
    msigdb_BP <- msigdbr(species = "Homo sapiens", category = "C5", subcategory = "BP")[, c("gs_name", "gene_symbol")]
    msigdb_MF <- msigdbr(species = "Homo sapiens", category = "C5", subcategory = "MF")[, c("gs_name", "gene_symbol")]

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
            axis.text.y = element_text(size = 14),
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
            axis.text.y = element_text(size = 14),
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
    # pdf(file = paste0("output/analysis/gsea_dotplots/run3/", iteration, "_", comp, "_dotplot.pdf"), width = 17, height = 15)
    # print(patchwork::wrap_plots(updotplot, downdotplot, ncol = 1))
    # dev.off()
    iteration <- iteration + 1
}

# Check for overlap between leading edge genes and rownames of the regulatory network
leading_edge_genes <- top_pathways_list$BP_blood_platelet_vs_lung_2pneumo$up[1, ]$leadingEdge[[1]]
# overlap_genes <- intersect(leading_edge_genes, colnames(lung_output$type_ii_pneumocyte$regNet))

# Extract the columns of the regulatory network that match the overlapping genes
overlap_columns <- lung_output$type_ii_pneumocyte$regNet[, leading_edge_genes, drop = FALSE]

# Save the top 10 highest values (and corresponding rows) for each column
top_edges <- apply(overlap_columns, 2, function(col) {
    sorted_indices <- order(col, decreasing = TRUE)[1:5]
    data.frame(tf = rownames(overlap_columns)[sorted_indices], edge = col[sorted_indices])
})

# Count the occurrences of each tf in the top_edges list
tf_counts <- sort(table(unlist(lapply(top_edges, function(df) df$tf))), decreasing = TRUE)
top_tfs <- names(tf_counts)[1:7]

lung_subset <- subset(lung_object, idents = c("type_ii_pneumocyte")) # load lung expression object first
blood_subset <- subset(blood_object, idents = c("platelet")) # load blood expression object first

# Extract the rows from the object and sum up expression
# overlap_sparse_matrix <- Seurat::AggregateExpression(lung_subset)[[1]][top_tfs, , drop = FALSE]
as.data.frame(rowSums(lung_subset@assays$RNA@counts[top_tfs, ]))
summary(rowSums(lung_subset@assays$RNA@counts))

# Print the resulting sparse matrix
# print(overlap_sparse_matrix)






#### 5. Run analysis on expression data ####

blood_object <- readRDS("output/preprocessing/correct_input_run3/blood_filter.rds")
Idents(blood_object) <- "cell_type"
lung_object <- readRDS("output/preprocessing/correct_input_run3/lung_filter.rds")
Idents(lung_object) <- "cell_type"

# blood_object[["RNA"]]@data <- sweep(blood_object[["RNA"]]@counts, 2, colSums(blood_object[["RNA"]]@counts), "/") * 1000000
# lung_object[["RNA"]]@data <- sweep(lung_object[["RNA"]]@counts, 2, colSums(lung_object[["RNA"]]@counts), "/") * 1000000

##### Calculate sum of expression for each gene in each cell type #####

calculate_expression_sum <- function(object) {
    sums <- list()
    for (type in levels(Idents(object))) {
        subset_object <- subset(object, idents = type)
        sum_exp <- rowSums(subset_object@assays$RNA@counts) # Take summed expression of gene per celltype
        sums[[type]] <- sum_exp
    }
    sums <- do.call(cbind, sums) # Combine list elements into a data frame structure
    return(sums)
}

blood_sum_exp <- calculate_expression_sum(blood_object)
lung_sum_exp <- calculate_expression_sum(lung_object)

colnames(blood_sum_exp) <- paste0("blood_", colnames(blood_sum_exp))
colnames(lung_sum_exp) <- paste0("lung_", colnames(lung_sum_exp))
combined_exp_matrix <- merge(blood_sum_exp, lung_sum_exp, by = "row.names", all = TRUE)
colnames(combined_exp_matrix)[1] <- "gene"
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
    gex_plot[[comp$name]] <- plot_linear_regression(data, comp$x, comp$y, comp$name, result$residuals, result$fit, "summed expression")
}

# Save plots of the linear regression models
ggsave("output/analysis/linear_regression/gex/immune_cells_between_tissues.pdf", patchwork::wrap_plots(A = gex_plot[[1]], B = gex_plot[[2]], C = gex_plot[[3]], design = "AABB\n#CC#"), width = 12, height = 6)
ggsave("output/analysis/linear_regression/gex/immune_cells_lung.pdf", patchwork::wrap_plots(A = gex_plot[[4]], B = gex_plot[[5]], C = gex_plot[[6]], design = "AABB\n#CC#"), width = 12, height = 6)
ggsave("output/analysis/linear_regression/gex/tissue_specific_linear.pdf", gex_plot[[7]], width = 12, height = 6)
ggsave("output/analysis/linear_regression/gex/blood_monocyte_vs_lung_2pneumo.pdf", gex_plot[[8]], width = 12, height = 6)

# Run GSEA on the expression data
gsea_results_gex <- run_gsea(gex_residuals)

# Generate plots for all comparisons in gsea_results and save to a list
gsea_gex_plots <- list()
gex_top_pathways_list <- list()
iteration <- 1
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
