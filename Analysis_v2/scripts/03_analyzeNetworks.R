#!/usr/bin/env Rscript


# Install required packages
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
if (!requireNamespace("org.Hs.eg.db", quietly = TRUE)) BiocManager::install("org.Hs.eg.db")
if (!requireNamespace("biomaRt", quietly = TRUE)) BiocManager::install("biomaRt")
if (!requireNamespace("msigdbr", quietly = TRUE)) BiocManager::install("msigdbr")

if (!requireNamespace("clusterProfiler", quietly = TRUE)) BiocManager::install("clusterProfiler")
if (!requireNamespace("AnnotationDbi", quietly = TRUE)) BiocManager::install("AnnotationDbi")

if (!requireNamespace("ggplot2", quietly = TRUE)) install.packages("ggplot2")

# Load required packages
library(ggplot2)
library(org.Hs.eg.db)
library(biomaRt)
library(msigdbr)
library(clusterProfiler) # gseGO
library(AnnotationDbi)

# Set working directory
setwd("/div/pythagoras/u1/siepv/siep/Analysis_v2/output/networks")

#### 1. Load and format data ####
load("final/bloodScorpionOutput.Rdata")
load("final/lungScorpionOutput.Rdata")
load("final/fatScorpionOutput.Rdata")
load("final/kidneyScorpionOutput.Rdata")


#### 2. Calculate indegrees & outdegrees ####
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



blood_outdegrees <- data.frame(tf = rownames(blood_output[[1]]$regNet))
for (celltype in names(blood_output)) {
    regnet <- blood_output[[celltype]]$regNet
    blood_outdegrees[[paste0("blood_", celltype)]] <- rowSums(regnet)
}


# Combine blood and lung indegrees and outdegrees into one dataframe
combined_indegrees <- Reduce(function(x, y) merge(x, y, by = "gene", all = TRUE), list(blood_indegrees, lung_indegrees, fat_indegrees, kidney_indegrees))

# Calculate correlation matrix
cor_indegree_matrix <- cor(combined_indegrees[, -1])

# Calculate RMSE??
rmse_indegree <- sqrt(mean((combined_indegrees$blood_cd4_positive_t_cell - combined_indegrees$lung_type_ii_pneumocyte)^2))



#### 3. Plot indegrees between conditions ####
# ggplot(combined_indegrees, aes(x = blood_cd4_positive_t_cell, y = blood_cd8_positive_t_cell)) +
#     geom_point() +
#     labs(
#         title = "Indegree Comparison between CD4+ and CD8+ T Cells in Blood",
#         x = "CD4+ T Cell Indegree",
#         y = "CD8+ T Cell Indegree"
#     ) +
#     geom_smooth(method = lm, formula = y ~ x, color = "red", fill = "#69b3a2", se = TRUE) +
#     theme_minimal()


#### 4. Calculate actual linear regression and residuals ####
x <- combined_indegrees$blood_cd4_positive_t_cell
y <- combined_indegrees$lung_type_ii_pneumocyte
fit <- lm(y ~ x, data = combined_indegrees)

# Extract residuals and rank them
res_df <- residuals(fit)
names(res_df) <- combined_indegrees$gene
ranked_df <- res_df[order(-res_df)]

# Plot the linear model and residuals
ggplot(combined_indegrees, aes(x = blood_platelet, y = lung_type_ii_pneumocyte)) +
    geom_point() +
    geom_smooth(method = "lm", formula = y ~ x, color = "red", fill = "#69b3a2", se = TRUE) +
    labs(
        title = "Linear Regression of CD4+ T Cell (blood) and type II pneumocytes Indegrees",
        x = "CD4+ T Cell Indegree",
        y = "Type II Pneumocyte Indegree"
    ) +
    theme_minimal()

# Plot residuals
residuals_df <- data.frame(gene = names(res_df), residuals = res_df)
ggplot(residuals_df, aes(x = seq_along(residuals), y = residuals)) +
    geom_point() +
    geom_hline(yintercept = 0, color = "red") +
    labs(
        title = "Residuals of Linear Regression",
        x = "ID",
        y = "Residuals"
    ) +
    theme_minimal()



############## Optie 1 #############

# Convert ENSEMBL IDs to ENTREZ IDs ### ensembl IDs worden niet herkend door gseGO?? dus eerst omzetten naar entrez
gene_conversion <- clusterProfiler::bitr(names(ranked_df), fromType = "ENSEMBL", toType = "ENTREZID", OrgDb = "org.Hs.eg.db")
ranked_list <- ranked_list[gene_conversion$ENSEMBL]
names(ranked_list) <- gene_conversion$ENTREZID

#### 5. Perform gene set enrichment analysis ####
gse_bp <- clusterProfiler::gseGO(ranked_df,
    ont = "BP",
    keyType = "ENSEMBL",
    OrgDb = "org.Hs.eg.db",
    pvalueCutoff = 0.5 # ,
    # eps = 1e-300
)

hist(gse_bp@result$p.adjust)

as.data.frame(gse_bp)


gse_mf <- clusterProfiler::gseGO(ranked_df,
    ont = "MF",
    keyType = "ENSEMBL",
    OrgDb = "org.Hs.eg.db",
    pvalueCutoff = 0.05 # ,
    # eps = 1e-300
)

as.data.frame(gse_mf)

# check hoeveel ID's er in de database overeenkomen
sum(names(ranked_df) %in% keys(org.Hs.eg.db, keytype = "ENSEMBL"))




########### Optie 2 #############

msigdb_BP <- msigdbr(species = "Homo sapiens", category = "C5", subcategory = "BP")
msigdb_BP <- msigdb_BP[, c("gs_name", "ensembl_gene")]
# does not work, term2gene needs a df
# msigdb_BP <- split(msigdb_BP$ensembl_gene, msigdb_BP$gs_name)

gsea_BP <- clusterProfiler::GSEA(ranked_df, TERM2GENE = msigdb_BP, pvalueCutoff = 1)

sum(names(ranked_df) %in% msigdb_BP$ensembl_gene)

converted <- clusterProfiler::bitr(names(ranked_df), fromType = "ENSEMBL", toType = "ENTREZID", OrgDb = "org.Hs.eg.db")





########### optie 3?? ##############

# Connect to the ENSEMBL database
ensembl <- biomaRt::useMart("ensembl", dataset = "hsapiens_gene_ensembl")

# Try to get updated mappings
converted_ids <- getBM(
    attributes = c("ensembl_gene_id", "ensembl_gene_id_version"),
    filters = "ensembl_gene_id",
    values = names(ranked_df),
    mart = ensembl
)
converted_ids <- converted_ids[!is.na(converted_ids$ensembl_gene_id_version), ]
ranked_df <- ranked_df[converted_ids$ensembl_gene_id]
names(ranked_df) <- converted_ids$ensembl_gene_id_version

gsea_BP <- clusterProfiler::GSEA(ranked_df, TERM2GENE = msigdb_BP, pvalueCutoff = 0.05)
