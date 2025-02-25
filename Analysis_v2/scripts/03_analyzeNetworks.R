#!/usr/bin/env Rscript


# Install required packages
if (!requireNamespace("ggplot2", quietly = TRUE)) install.packages("ggplot2")
if (!requireNamespace("Matrix", quietly = TRUE)) install.packages("Matrix")

if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
if (!requireNamespace("msigdbr", quietly = TRUE)) BiocManager::install("msigdbr")
if (!requireNamespace("fgsea", quietly = TRUE)) BiocManager::install("fgsea")

# Load required packages
library(ggplot2)
library(msigdbr)
library(fgsea)
library(Matrix)

# Set working directory
setwd("/div/pythagoras/u1/siepv/siep/Analysis_v2/output/networks")

#### 1. Load and format data ####
load("final/bloodScorpionOutput.Rdata")
load("final/lungScorpionOutput.Rdata")
load("final/fatScorpionOutput.Rdata")
load("final/kidneyScorpionOutput.Rdata")
load("final/liverScorpionOutput.Rdata")


#### 2. Calculate indegrees & outdegrees ####       #rewrite to function!!
blood_indegrees <- data.frame(gene = colnames(blood_output[[1]]$regNet))
for (celltype in names(blood_output)) {
    regnet <- blood_output[[celltype]]$regNet
    blood_indegrees[[paste0("blood_", celltype)]] <- colSums(regnet)
}
lung_indegrees <- data.frame(gene = colnames(lung_output[[1]]$regNet))
for (celltype in names(lung_output)) {
    regnet <- lung_output[[celltype]]$regNet
    lung_indegrees[[paste0("lung_", celltype)]] <- colSums(regnet)
}
fat_indegrees <- data.frame(gene = colnames(fat_output[[1]]$regNet))
for (celltype in names(fat_output)) {
    regnet <- fat_output[[celltype]]$regNet
    fat_indegrees[[paste0("fat_", celltype)]] <- colSums(regnet)
}
kidney_indegrees <- data.frame(gene = colnames(kidney_output[[1]]$regNet))
for (celltype in names(kidney_output)) {
    regnet <- kidney_output[[celltype]]$regNet
    kidney_indegrees[[paste0("kidney_", celltype)]] <- colSums(regnet)
}
liver_indegrees <- data.frame(gene = colnames(liver_output[[1]]$regNet))
for (celltype in names(liver_output)) {
    regnet <- liver_output[[celltype]]$regNet
    liver_indegrees[[paste0("liver_", celltype)]] <- colSums(regnet)
}

# Combine blood and lung indegrees and outdegrees into one dataframe
combined_indegrees <- Reduce(function(x, y) merge(x, y, by = "gene", all = TRUE), list(blood_indegrees, lung_indegrees, fat_indegrees, kidney_indegrees, liver_indegrees))
# combined_indegrees <- na.omit(combined_indegrees)

# Calculate correlation matrix
cor_indegree_matrix <- cor(combined_indegrees[, -1])


#### 3. Calculate linear regression and residuals between conditions ####

# Define comparisons
comparisons <- list(
    list(name = "blood_platelet_vs_lung_pneumo", x = "blood_platelet", y = "lung_type_ii_pneumocyte"),
    list(name = "blood_cd4_vs_cd8", x = "blood_cd4_positive_t_cell", y = "blood_cd8_positive_t_cell"),
    list(name = "blood_cd4_vs_monocyte", x = "blood_cd4_positive_t_cell", y = "blood_monocyte"),
    list(name = "lung_cd4_vs_cd8", x = "lung_cd4_positive_t_cell", y = "lung_cd8_positive_t_cell"),
    list(name = "lung_cd4_vs_monocyte", x = "lung_cd4_positive_t_cell", y = "lung_monocyte"),
    list(name = "cd4_blood_vs_lung", x = "blood_cd4_positive_t_cell", y = "lung_cd4_positive_t_cell"),
    list(name = "cd8_blood_vs_lung", x = "blood_cd8_positive_t_cell", y = "lung_cd8_positive_t_cell")
)

# Function to fit linear model and extract + rank residuals descending
calculate_residuals <- function(data, x_con, y_con) {
    fit <- lm(data[[y_con]] ~ data[[x_con]], data = data) # Fit linear model
    res <- residuals(fit)
    names(res) <- data$gene
    return(res[order(-res)]) # Rank residuals
}

residuals_list <- list()
for (comp in comparisons) {
    cat("Calculating residuals for:", comp$name, "\n")
    residuals_list[[comp$name]] <- calculate_residuals(combined_indegrees, comp$x, comp$y)
}


#### 4. Run GSEA on residuals ####

# Load gene sets for GO:BP and GO:MF using gene symbols
msigdb_BP <- msigdbr(species = "Homo sapiens", category = "C5", subcategory = "BP")[, c("gs_name", "gene_symbol")]
msigdb_MF <- msigdbr(species = "Homo sapiens", category = "C5", subcategory = "MF")[, c("gs_name", "gene_symbol")]

# Convert to fgsea-compatible format (list of gene sets)
msigdb_BP_list <- split(msigdb_BP$gene_symbol, msigdb_BP$gs_name)
msigdb_MF_list <- split(msigdb_MF$gene_symbol, msigdb_MF$gs_name)


gsea_results <- list()

for (comp in names(residuals_list[1:2])) {
    cat("Running GSEA for:", comp, "\n")

    ranked_genes <- residuals_list[[comp]] # Get ranked gene list

    # Run fgsea for GO:BP
    gsea_results[[paste0(comp, "_BP")]] <- fgsea(
        pathways = msigdb_BP_list,
        stats = ranked_genes,
        minSize = 10,
        maxSize = 500,
        nproc = 1
    )

    # Run fgsea for GO:MF
    gsea_results[[paste0(comp, "_MF")]] <- fgsea(
        pathways = msigdb_MF_list,
        stats = ranked_genes,
        minSize = 10,
        maxSize = 500,
        nperm = 10000
    )
}





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



top_genes <- names(residuals_list$blood_platelet_vs_lung_pneumo)[1:1000] # Take the top 500 most deviated genes
neg_genes <- names(rev(residuals_list$blood_platelet_vs_lung_pneumo))[1:1000]

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
    values = missing_top_genes, # or missing_neg_genes
    mart = ensembl
)

# Check how many are lncRNAs
lncRNAs <- subset(gene_annotations, gene_biotype == "lncRNA")
prcoding <- subset(gene_annotations, gene_biotype == "protein_coding")
head(prcoding)
