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
library(stringr)
library(Seurat)
# library(biomaRt)

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

# Calculate correlation matrix
cor_indegree_matrix <- cor(combined_indegrees[, -1])

# Convert the correlation matrix to a long format for ggplot
cor_long <- as.data.frame(as.table(cor_indegree_matrix))

setwd("/div/pythagoras/u1/siepv/siep/Analysis_v2/output/analysis")

# Plot the one-sided correlation matrix with values in the tiles
ggsave("correlation/correlation_matrix_indegrees.pdf", ggplot(cor_long, aes(Var1, Var2, fill = Freq)) +
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

# Function to plot linear regression with ggplot and label extreme residuals
plot_linear_regression <- function(data, x_con, y_con, comp_name, residuals) {
    top_gene <- names(residuals)[1] # Gene with highest residual
    bottom_gene <- names(residuals)[length(residuals)] # Gene with lowest residual

    ggplot(data, aes(x = .data[[x_con]], y = .data[[y_con]])) +
        geom_point(alpha = 0.5) +
        geom_point(
            data = data %>% filter(gene %in% c(top_gene, bottom_gene)),
            aes(color = gene),
            size = 3
        ) +
        scale_color_manual(values = setNames(c("red", "blue"), c(top_gene, bottom_gene))) +
        geom_smooth(method = "lm", se = FALSE, color = "red") +
        geom_text(
            data = data %>% filter(gene == top_gene),
            aes(label = gene, color = gene),
            vjust = -1,
            size = 5,
            fontface = "bold"
        ) +
        geom_text(
            data = data %>% filter(gene == bottom_gene),
            aes(label = gene, color = gene),
            vjust = 1.5,
            size = 5,
            fontface = "bold"
        ) +
        labs(
            title = paste("Linear Regression:", comp_name),
            x = x_con,
            y = y_con
        ) +
        theme_minimal() +
        theme(legend.position = "none")
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
    plot_list[[comp$name]] <- plot_linear_regression(combined_indegrees, comp$x, comp$y, comp$name, result$residuals)
}

# Save plots of the linear regression models
ggsave("linear_regression/immune_cells_between_tissues.pdf", patchwork::wrap_plots(A = plot_list[[1]], B = plot_list[[2]], C = plot_list[[3]], design = "AABB\n#CC#"), width = 12, height = 6)
ggsave("linear_regression/immune_cells_lung.pdf", patchwork::wrap_plots(A = plot_list[[4]], B = plot_list[[5]], C = plot_list[[6]], design = "AABB\n#CC#"), width = 12, height = 6)
ggsave("linear_regression/tissue_specific_linear.pdf", plot_list[[7]], width = 12, height = 6)



#### 4. Run GSEA on ranked residuals ####

# Load gene sets for GO:BP and GO:MF using gene symbols
msigdb_BP <- msigdbr(species = "Homo sapiens", category = "C5", subcategory = "BP")[, c("gs_name", "gene_symbol", "ensembl_gene")]
msigdb_MF <- msigdbr(species = "Homo sapiens", category = "C5", subcategory = "MF")[, c("gs_name", "gene_symbol", "ensembl_gene")]

# Convert to fgsea-compatible format (list of gene sets)
msigdb_BP <- split(msigdb_BP$gene_symbol, msigdb_BP$gs_name)
msigdb_MF <- split(msigdb_MF$gene_symbol, msigdb_MF$gs_name)

msigdb_BP_combined <- split(c(msigdb_BP$gene_symbol, msigdb_BP$ensembl_gene), msigdb_BP$gs_name)


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

# Function to get top up and down pathways
get_top_pathways <- function(gsea_result) {
    sig_pathways <- gsea_result %>%
        filter(padj < 0.05)

    sig_pathways <- sig_pathways %>%
        mutate(
            pathway = gsub("_", " ", as.character(pathway)),
            pathway = str_wrap(pathway, width = 50, whitespace_only = FALSE),
            pathway = gsub(" ", "_", pathway)
        )

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

    return(list(up = topPathwaysUp, down = topPathwaysDown))
}

# Function to plot top pathways
plot_top_pathways <- function(topPathways, comparison_name) {
    p1 <- ggplot(topPathways$up, aes(x = NES, y = reorder(pathway, NES), size = size, fill = NES)) +
        geom_point(shape = 21) +
        scale_size_area(max_size = 10) +
        scale_fill_gradient(low = "white", high = "red") +
        theme(
            axis.text.y = element_text(size = 14),
            axis.title.y = element_blank(),
            plot.title = element_text(hjust = 1, face = "bold") # Make title bold
        ) +
        xlim(min(topPathways$up$NES), max(topPathways$up$NES)) +
        xlab("NES") +
        ggtitle(paste("Top 10 positively enriched pathways in", comparison_name))

    p2 <- ggplot(topPathways$down, aes(x = NES, y = pathway, size = size, fill = NES)) +
        geom_point(shape = 21) +
        scale_size_area(max_size = 10) +
        scale_fill_gradient(low = "blue", high = "white") +
        theme(
            axis.text.y = element_text(size = 14),
            axis.title.y = element_blank(),
            plot.title = element_text(hjust = 1, face = "bold") # Make title bold
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
for (comp in names(gsea_results)) {
    # Generate plots
    cat("Plotting GSEA results for:", comp, "\n")
    topPathways <- get_top_pathways(gsea_results[[comp]])
    gsea_plots[[comp]] <- plot_top_pathways(topPathways, comp)
    top_pathways_list[[comp]] <- topPathways

    # Save plots
    updotplot <- gsea_plots[[comp]]$up
    downdotplot <- gsea_plots[[comp]]$down
    pdf(file = paste0("gsea_dotplots/", iteration, "_", comp, "_dotplot.pdf"), width = 17, height = 15)
    print(patchwork::wrap_plots(updotplot, downdotplot, ncol = 1))
    dev.off()
    iteration <- iteration + 1
}



# pdf(file = "13_enrichmentPlot_last.pdf", width = 20, height = 15)
# plotEnrichment(
#     msigdb_BP[[gsea_results[[13]][order(NES)][2]$pathway]],
#     ranked_genes
# ) +
#     labs(title = gsea_results[[13]][order(NES)][2]$pathway)
# dev.off()



#### 5. Run GSEA on expression data ####

setwd("/div/pythagoras/u1/siepv/siep/Analysis_v2")

blood_object <- readRDS("output/preprocessing/blood_filter.rds")
Idents(blood_object) <- "cell_type"
lung_object <- readRDS("output/preprocessing/lung_filter.rds")
Idents(lung_object) <- "cell_type"

##### Calculate sum of expression for each gene in each cell type #####
blood_sum_exp <- Seurat::AggregateExpression(blood_object)[[1]]
lung_sum_exp <- Seurat::AggregateExpression(lung_object)[[1]]

colnames(blood_sum_exp) <- paste0("blood_", colnames(blood_sum_exp))
colnames(lung_sum_exp) <- paste0("lung_", colnames(lung_sum_exp))
combined_exp_matrix <- merge(blood_sum_exp, lung_sum_exp, by = "row.names", all = TRUE)
colnames(combined_exp_matrix)[1] <- "gene"
cor_exp_matrix <- cor(combined_exp_matrix[, -1])

##### Calculate residuals for each comparison #####
gex_residuals <- list()
gex_fit <- list()
gex_plot <- list()
for (comp in comparisons) {
    cat("Calculating gex residuals for:", comp$name, "\n")
    result <- calculate_residuals(combined_exp_matrix, comp$x, comp$y)
    gex_residuals[[comp$name]] <- result$residuals
    gex_fit[[comp$name]] <- result$fit
    gex_plot[[comp$name]] <- plot_linear_regression(combined_exp_matrix, comp$x, comp$y, comp$name, result$residuals)
}

table(stringr::str_extract(names(gex_residuals[[5]]), "ENSG[0-9]+.[0-9]+"))

table(grepl("ENSG[0-9]+.[0-9]+", names(gex_residuals[[5]])))

table(grepl("^ENSG[0-9]+\\.[0-9]+$", names(gex_residuals[[5]])))











#### GSEA testing stuff ####

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
