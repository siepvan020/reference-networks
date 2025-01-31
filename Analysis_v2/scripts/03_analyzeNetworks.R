#!/usr/bin/env Rscript


# Install required packages
if (!requireNamespace("Seurat", quietly = TRUE)) install.packages("Seurat")

# Load required packages
library(Seurat)
library(ggplot2)

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
ggplot(combined_indegrees, aes(x = blood_cd4_positive_t_cell, y = blood_cd8_positive_t_cell)) +
    geom_point() +
    labs(
        title = "Indegree Comparison between CD4+ and CD8+ T Cells in Blood",
        x = "CD4+ T Cell Indegree",
        y = "CD8+ T Cell Indegree"
    ) +
    geom_smooth(method = lm, formula = y ~ x, color = "red", fill = "#69b3a2", se = TRUE) +
    theme_minimal()


#### 4. Calculate actual linear regression ####
# x <-
# y <-
# lm(y ~ x, )
