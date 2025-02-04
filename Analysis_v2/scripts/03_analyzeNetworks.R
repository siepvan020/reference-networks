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


#### 2. Calculate indegrees ####
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

# Combine blood and lung indegrees into one dataframe
combined_indegrees <- merge(blood_indegrees, lung_indegrees, by = "gene")


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
x <- blood_indegrees$blood_cd4_positive_t_cell
y <- blood_indegrees$blood_cd8_positive_t_cell
fit <- lm(y ~ x, data = blood_indegrees)

# Extract residuals and rank them
res_df <- residuals(fit)
names(res_df) <- combined_indegrees$gene
ranked_df <- res_df[order(-res_df)]


############## Optie 1 #############

# Convert ENSEMBL IDs to ENTREZ IDs ### ensembl IDs worden niet herkend door gseGO?? dus eerst omzetten naar entrez
# gene_conversion <- clusterProfiler::bitr(names(ranked_df), fromType = "ENSEMBL", toType = "ENTREZID", OrgDb = "org.Hs.eg.db")
# ranked_list <- ranked_list[gene_conversion$ENSEMBL]
# names(ranked_list) <- gene_conversion$ENTREZID

# ggplot(fit_df, aes(x = fitted, y = residuals)) +
#     geom_point() +
#     geom_hline(yintercept = 0, color = "red") +
#     labs(
#         title = "Residuals vs Fitted Values",
#         x = "Fitted Values",
#         y = "Residuals"
#     ) +
#     theme_minimal()

#### 5. Perform gene set enrichment analysis ####
gse_bp <- clusterProfiler::gseGO(ranked_df,
    ont = "BP",
    keyType = "ENSEMBL",
    OrgDb = "org.Hs.eg.db",
    eps = 1e-300
)

as.data.frame(gse_bp)

gse_mf <- clusterProfiler::gseGO(ranked_df,
    ont = "MF",
    keyType = "ENSEMBL",
    OrgDb = "org.Hs.eg.db" # ,
    # eps = 1e-300
)

# check hoeveel ID's er in de database overeenkomen
sum(names(ranked_df) %in% keys(org.Hs.eg.db, keytype = "ENSEMBL"))



########### Optie 2 #############

msigdb_MF <- msigdbr(species = "Homo sapiens", category = "C5", subcategory = "MF")
msigdb_MF <- msigdb_BP[, c("gs_name", "ensembl_gene")]
# does not work, term2gene needs a df
# msigdb_BP <- split(msigdb_BP$ensembl_gene, msigdb_BP$gs_name)

gsea_MF <- clusterProfiler::GSEA(ranked_df, TERM2GENE = msigdb_MF, pvalueCutoff = 0.05)

sum(names(ranked_df) %in% msigdb_BP$ensembl_gene)

converted <- clusterProfiler::bitr(names(ranked_df), fromType = "ENSEMBL", toType = "ENTREZID", OrgDb = "org.Hs.eg.db")





########### optie 3?? ##############

# Connect to the ENSEMBL database
ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

# Check if your IDs exist in the latest ENSEMBL database
valid_genes <- getBM(
    attributes = c("ensembl_gene_id"),
    filters = "ensembl_gene_id",
    values = names(ranked_df),
    mart = ensembl
)

# Count how many of your IDs exist in ENSEMBL
sum(names(ranked_df) %in% valid_genes$ensembl_gene_id)
