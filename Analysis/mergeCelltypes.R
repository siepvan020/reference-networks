#!/usr/bin/env Rscript

# Set CRAN mirror
options(repos = c(CRAN = "https://cloud.r-project.org"))

# Load required packages
if (!requireNamespace("remotes", quietly = TRUE)) install.packages("remotes")
if (!requireNamespace("SeuratDisk", quietly = TRUE)) remotes::install_github("mojaveazure/seurat-disk")
if (!requireNamespace("SCORPION", quietly = TRUE)) install.packages("SCORPION")

library(SeuratDisk)
library(SCORPION)
library(Seurat)

# Set working directory
setwd("/div/pythagoras/u1/siepv/siep/Analysis")

# Function to load, set identities, and subset Seurat objects
process_seurat_object <- function(file_path, mapping) {
    setwd(dirname(file_path))
    seurat_object <- LoadH5Seurat(basename(file_path), assays = "decontXcounts", reductions = FALSE)

    Idents(seurat_object) <- "cell_ontology_class"
    Idents(seurat_object) <- factor(Idents(seurat_object), levels = names(mapping), labels = mapping)
    seurat_object <- subset(seurat_object, idents = unique(mapping))

    return(seurat_object)
}

# Define mappings for cell ontology classes to simplified cell types
lung_mapping <- c(
    "cd4-positive alpha-beta t cell" = "cd4_positive_t_cell",
    "cd4-positive, alpha-beta t cell" = "cd4_positive_t_cell",
    "cd8-positive alpha-beta t cell" = "cd8_positive_t_cell",
    "cd8-positive, alpha-beta t cell" = "cd8_positive_t_cell",
    "classical monocyte" = "monocyte",
    "intermediate monocyte" = "monocyte",
    "non-classical monocyte" = "monocyte",
    "type ii pneumocyte" = "type_ii_pneumocyte"
)

blood_mapping <- c(
    "cd4-positive, alpha-beta memory t cell" = "cd4_positive_t_cell",
    "cd4-positive, alpha-beta t cell" = "cd4_positive_t_cell",
    "naive thymus-derived cd4-positive, alpha-beta t cell" = "cd4_positive_t_cell",
    "cd8-positive, alpha-beta cytokine secreting effector t cell" = "cd8_positive_t_cell",
    "cd8-positive, alpha-beta t cell" = "cd8_positive_t_cell",
    "classical monocyte" = "monocyte",
    "monocyte" = "monocyte",
    "non-classical monocyte" = "monocyte",
    "platelet" = "platelet"
)

# Process lung and blood Seurat objects
print("LUNG")
lung_object <- process_seurat_object("/div/pythagoras/u1/siepv/siep/Analysis/data/h5seurat/TS_Lung.h5seurat", lung_mapping)
print("BLOOD")
blood_object <- process_seurat_object("/div/pythagoras/u1/siepv/siep/Analysis/data/h5seurat/TS_Blood.h5seurat", blood_mapping)

print("SAVE")
setwd("/div/pythagoras/u1/siepv/siep/Analysis/Rdata")
saveRDS(lung_object, "lung_object.rds")
saveRDS(blood_object, "blood_object.rds")


lung_object <- LoadH5Seurat("/div/pythagoras/u1/siepv/siep/Analysis/data/h5seurat/TS_Lung.h5seurat", misc = FALSE)
blood_object <- LoadH5Seurat("/div/pythagoras/u1/siepv/siep/Analysis/data/h5seurat/TS_Blood.h5seurat", misc = FALSE)

DimPlot(lung_object, reduction = "umap", group.by = "cell_ontology_class", label = TRUE)
DimPlot(blood_object, reduction = "umap", group.by = "cell_ontology_class", label = TRUE)
