#!/usr/bin/env Rscript


# Install required packages
if (!requireNamespace("Seurat", quietly = TRUE)) install.packages("Seurat")

# Load required packages
library(Seurat)
library(ggplot2)

# Set working directory
setwd("/div/pythagoras/u1/siepv/siep/Analysis_v2/output/networks")

#### Load and format data ####
load("final/bloodScorpionOutput.Rdata")
blood_output <- scorpion_output
load("final/lungScorpionOutput.Rdata")
lung_output <- scorpion_output

#### Add indegrees and outdegrees####         gives wrong/repetitive degree values, fix later
blood_dfs <- list()
for (celltype in names(blood_output)) {
    regnet <- blood_output[[celltype]]$regNet
    df <- as.data.frame(as.table(as(regnet, "matrix")))
    colnames(df) <- c("TF", "Target_Gene", "Weight")
    df$Indegree <- colSums(regnet)
    df$Outdegree <- rowSums(regnet)
    blood_dfs[[celltype]] <- df
}

lung_dfs <- list() #   gives wrong/repetitive degree values, fix later
for (celltype in names(lung_output)) {
    regnet <- lung_output[[celltype]]$regNet
    df <- as.data.frame(as.table(as(regnet, "matrix")))
    colnames(df) <- c("TF", "Target_Gene", "Weight")
    df$Indegree <- colSums(regnet)
    df$Outdegree <- rowSums(regnet)
    lung_dfs[[celltype]] <- df
}


# temp <- colSums(lung_output$monocyte$regNet)
