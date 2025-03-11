#!/usr/bin/env Rscript


# Install required packages
if (!requireNamespace("ggplot2", quietly = TRUE)) install.packages("ggplot2")
if (!requireNamespace("Matrix", quietly = TRUE)) install.packages("Matrix")

if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
if (!requireNamespace("msigdbr", quietly = TRUE)) BiocManager::install("msigdbr")
if (!requireNamespace("fgsea", quietly = TRUE)) BiocManager::install("fgsea")

# Load required packages
library(ggplot2)
library(patchwork)
library(msigdbr)
library(fgsea)
library(Matrix)
library(dplyr)
library(biomaRt)

# Set working directory
setwd("/div/pythagoras/u1/siepv/siep/Analysis_v2/output/networks")

#### 1. Load and format data ####
load("final/bloodScorpionOutput.Rdata")
load("final/lungScorpionOutput.Rdata")
load("final/fatScorpionOutput.Rdata")
load("final/kidneyScorpionOutput.Rdata")
load("final/liverScorpionOutput.Rdata")


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
fat_indegrees <- calculate_indegrees(fat_output, "fat")
kidney_indegrees <- calculate_indegrees(kidney_output, "kidney")
liver_indegrees <- calculate_indegrees(liver_output, "liver")

# Combine blood and lung indegrees and outdegrees into one dataframe
combined_indegrees <- Reduce(function(x, y) merge(x, y, by = "gene", all = TRUE), list(blood_indegrees, lung_indegrees)) # , fat_indegrees, kidney_indegrees, liver_indegrees))
# combined_indegrees <- na.omit(combined_indegrees)

# Calculate correlation matrix
cor_indegree_matrix <- cor(combined_indegrees[, -1])

# Convert the correlation matrix to a long format for ggplot
cor_long <- as.data.frame(as.table(cor_indegree_matrix))

setwd("/div/pythagoras/u1/siepv/siep/Analysis_v2/output/analysis")

# Plot the one-sided correlation matrix with values in the tiles
ggsave("correlation_matrix_indegrees.pdf", ggplot(cor_long, aes(Var1, Var2, fill = Freq)) +
    geom_tile(na.rm = TRUE) +
    geom_text(aes(label = round(Freq, 3)), color = "black", size = 3, na.rm = TRUE) +
    scale_fill_gradient(low = "white", high = "steelblue", name = "Pearson's r value") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
    labs(x = "Cell Types", y = "Cell Types", title = "Correlation Matrix of Indegrees"),
width = 12, height = 6
)

#### 3. Calculate linear regression and residuals between conditions ####

# Define comparisons
comparisons <- list(
    list(name = "cd4_blood_vs_lung", x = "blood_cd4_positive_t_cell", y = "lung_cd4_positive_t_cell"),
    list(name = "cd8_blood_vs_lung", x = "blood_cd8_positive_t_cell", y = "lung_cd8_positive_t_cell"),
    list(name = "monocyte_blood_vs_lung", x = "blood_monocyte", y = "lung_monocyte"),
    list(name = "lung_cd4_vs_cd8", x = "lung_cd4_positive_t_cell", y = "lung_cd8_positive_t_cell"),
    list(name = "lung_cd4_vs_monocyte", x = "lung_cd4_positive_t_cell", y = "lung_monocyte"),
    list(name = "lung_cd8_vs_monocyte", x = "lung_cd8_positive_t_cell", y = "lung_monocyte"),
    list(name = "blood_platelet_vs_lung_2pneumo", x = "blood_platelet", y = "lung_type_ii_pneumocyte")
)

# Function to fit linear model and extract + rank residuals descending
calculate_residuals <- function(data, x_con, y_con) {
    fit <- lm(data[[y_con]] ~ data[[x_con]], data = data) # Fit linear model
    res <- residuals(fit)
    names(res) <- data$gene
    return(list(fit = fit, residuals = res[order(-res)])) # Return fit and ranked residuals
}

# Function to plot linear regression with ggplot
plot_linear_regression <- function(data, x_con, y_con, comp_name) {
    ggplot(data, aes(x = x_con, y = y_con)) +
        geom_point(alpha = 0.5) +
        geom_smooth(method = "lm", se = FALSE, color = "red") +
        labs(
            title = paste("Linear Regression:", comp_name),
            x = x_con,
            y = y_con
        ) +
        theme_minimal()
}

# Calculate residuals for each comparison
residuals_list <- list()
fit_list <- list()
plot_list <- list()
for (comp in comparisons) {
    cat("Calculating residuals for:", comp$name, "\n")
    result <- calculate_residuals(combined_indegrees, comp$x, comp$y)
    residuals_list[[comp$name]] <- result$residuals
    fit_list[[comp$name]] <- result$fit
    plot_list[[comp$name]] <- plot_linear_regression(combined_indegrees, comp$x, comp$y, comp$name)
}

# Save plots of the linear regression models
ggsave("immune_cells_between_tissues.pdf", patchwork::wrap_plots(A = plot_list[[1]], B = plot_list[[2]], C = plot_list[[3]], design = "AABB\n#CC#"), width = 12, height = 6)
ggsave("immune_cells_lung.pdf", patchwork::wrap_plots(A = plot_list[[4]], B = plot_list[[5]], C = plot_list[[6]], design = "AABB\n#CC#"), width = 12, height = 6)
ggsave("tissue_specific_linear.pdf", plot_list[[7]], width = 12, height = 6)



#### 4. Run GSEA on ranked residuals ####

# Load gene sets for GO:BP and GO:MF using gene symbols
msigdb_BP <- msigdbr(species = "Homo sapiens", category = "C5", subcategory = "BP")[, c("gs_name", "gene_symbol")]
msigdb_MF <- msigdbr(species = "Homo sapiens", category = "C5", subcategory = "MF")[, c("gs_name", "gene_symbol")]

# Convert to fgsea-compatible format (list of gene sets)
msigdb_BP <- split(msigdb_BP$gene_symbol, msigdb_BP$gs_name)
msigdb_MF <- split(msigdb_MF$gene_symbol, msigdb_MF$gs_name)


gsea_results <- list()

options(warn = 1)
for (comp in names(residuals_list)) {
    cat("Running GSEA for:", comp, "\n")

    ranked_genes <- residuals_list[[comp]] # Get ranked gene list

    # Run fgsea for GO:BP
    gsea_results[[paste0("BP_", comp)]] <- fgsea(
        pathways = msigdb_BP,
        stats = ranked_genes,
        minSize = 10,
        maxSize = 500,
        nproc = 1
    )

    # Run fgsea for GO:MF
    gsea_results[[paste0("MF_", comp)]] <- fgsea(
        pathways = msigdb_MF,
        stats = ranked_genes,
        minSize = 10,
        maxSize = 500,
        nproc = 1
    )
}


# Function to plot top pathways
plot_top_pathways <- function(gsea_result, comparison_name) {
    sig_pathways <- gsea_result %>%
        filter(padj < 0.05)

    topPathwaysUp <- head(sig_pathways[NES > 0], 10) %>%
        arrange(NES)
    topPathwaysUp$pathway <- factor(topPathwaysUp$pathway, levels = topPathwaysUp$pathway)

    topPathwaysDown <- head(sig_pathways[NES < 0], 10) %>%
        arrange(NES)
    topPathwaysDown$pathway <- factor(topPathwaysDown$pathway, levels = topPathwaysDown$pathway)

    topPathways <- bind_rows(topPathwaysUp, topPathwaysDown) %>%
        arrange(NES)
    topPathways$pathway <- factor(topPathways$pathway, levels = topPathways$pathway)

    p1 <- ggplot(topPathwaysUp, aes(x = NES, y = pathway, size = size, fill = NES)) +
        geom_point(shape = 21) +
        scale_size_area(max_size = 10) +
        scale_fill_gradient(low = "white", high = "red") +
        theme(axis.text.y = element_text(size = 11), axis.title.y = element_blank()) +
        xlim(min(topPathwaysUp$NES), max(topPathwaysUp$NES)) +
        xlab("NES") +
        ggtitle(paste("Top 10 enriched pathways in", comparison_name, "upregulated"))

    p2 <- ggplot(topPathwaysDown, aes(x = NES, y = pathway, size = size, fill = NES)) +
        geom_point(shape = 21) +
        scale_size_area(max_size = 10) +
        scale_fill_gradient(low = "blue", high = "white") +
        theme(axis.text.y = element_text(size = 11), axis.title.y = element_blank()) +
        xlim(min(topPathwaysDown$NES), max(topPathwaysDown$NES)) +
        xlab("NES") +
        ggtitle(paste("Top 10 enriched pathways in", comparison_name, "downregulated"))

    return(list(up = p1, down = p2))
}

# Generate plots for all comparisons in gsea_results and save to a list
gsea_plots <- list()
for (comp in names(gsea_results)) {
    cat("Plotting GSEA results for:", comp, "\n")
    gsea_plots[[comp]] <- plot_top_pathways(gsea_results[[comp]], comp)
}

dotplot <- gsea_plots[[7]]$up
print(dotplot)
filename <- paste0(names(gsea_plots)[7], "_gsea_plot.pdf")
ggsave(filename, dotplot, width = 12, height = 8)



# pdf(file = "7_gsea_top20pathways.pdf", width = 20, height = 15)
# plotGseaTable(msigdb_BP[topPathways], stats = ranked_genes, fgseaRes = sig_pathways, gseaParam = 0.5)
# dev.off()


# plotEnrichment(
#     msigdb_BP_list[[gsea_results[[1]][order(padj), ][1]$pathway]],
#     ranked_genes
# ) +
#     labs(title = gsea_results[[1]][order(padj), ][1]$pathway)




















##### plot gsea test

library(ggplot2)
library(dplyr)

gsea_results$blood_platelet_vs_lung_pneumo_BP %>%
    arrange(padj) %>%
    head(20) %>%
    ggplot(aes(reorder(pathway, NES), NES, fill = padj < 0.05)) +
    geom_col() +
    coord_flip() +
    scale_fill_manual(values = c("TRUE" = "red", "FALSE" = "gray")) +
    labs(
        x = "Pathway",
        y = "Normalized Enrichment Score (NES)",
        title = "Top 20 Enriched Pathways"
    ) +
    theme_minimal()


# Extract the first pathway name
first_pathway <- gsea_results[[1]]$pathway[1]

# Extract the ranked gene list for the first comparison
ranked_genes <- residuals_list[[names(residuals_list)[1]]]
ranked_genes <- sort(ranked_genes, decreasing = TRUE) # Ensure it's sorted

# Plot the enrichment for the first pathway
plotEnrichment(
    pathway = msigdb_BP_list[["GOBP_TRNA_METABOLIC_PROCESS"]], # Get gene set for first pathway
    stats = ranked_genes
) + labs(title = paste("Enrichment Plot:", first_pathway))



top_genes <- names(residuals_list[[7]])[1:1000] # Take the top 1000 most deviated genes
neg_genes <- names(rev(residuals_list[[7]]))[1:1000]

missing_top_genes <- setdiff(top_genes, msigdb_BP$gene_symbol)
missing_neg_genes <- setdiff(neg_genes, msigdb_BP$gene_symbol)

mf_missing_top_genes <- setdiff(top_genes, msigdb_MF$gene_symbol)
mf_missing_neg_genes <- setdiff(neg_genes, msigdb_MF$gene_symbol)

length(missing_top_genes) # How many top genes are missing?
length(missing_neg_genes) # How many bottom genes are missing?
head(missing_top_genes) # See which ones
head(missing_neg_genes) # See which ones


library(biomaRt)

# Connect to Ensembl
ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

# Get gene biotypes for missing genes
gene_annotations <- getBM(
    attributes = c("external_gene_name", "gene_biotype"),
    filters = "external_gene_name",
    values = names(head(residuals_list[[1]], 1000)), # or missing_neg_genes
    mart = ensembl
)

# Check how many are lncRNAs
lncRNAs <- subset(gene_annotations, gene_biotype == "lncRNA")
prcoding <- subset(gene_annotations, gene_biotype == "protein_coding")
head(prcoding)
