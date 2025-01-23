#!/usr/bin/env Rscript


progress_file <- "/div/pythagoras/u1/siepv/siep/Analysis_v2/v2_parallel/network_progress.log"

# Start fresh log or append
cat("Starting script at", Sys.time(), "\n", file = progress_file)

# Install required packages
if (!requireNamespace("Seurat", quietly = TRUE)) install.packages("Seurat")
if (!requireNamespace("SCORPION", quietly = TRUE)) install.packages("SCORPION")
if (!requireNamespace("doParallel", quietly = TRUE)) install.packages("doParallel")
if (!requireNamespace("Matrix", quietly = TRUE)) install.packages("Matrix")

# Load required packages
library(Seurat)
library(SCORPION)
library(doParallel)
library(Matrix)

# Set working directory
setwd("/div/pythagoras/u1/siepv/siep/Analysis_v2")

# Load Seurat objects
blood_object <- readRDS("output/celltype_filtering/blood_celltypes.rds")
lung_object <- readRDS("output/celltype_filtering/lung_celltypes.rds")
cat("- Loaded Seurat objects- \n", file = progress_file, append = TRUE)

# Read priors
ppi <- read.delim("data/priors/ppi_prior_2024.tsv", header = FALSE, sep = "\t")
tf <- read.delim("data/priors/motif_prior_names_2024.tsv", header = FALSE, sep = "\t")
cat("- Loaded priors - \n", file = progress_file, append = TRUE)

# Set up parallel processing
num_cores <- 4
cl <- makeCluster(num_cores)
registerDoParallel(cl)
clusterExport(cl, varlist = c("tf", "ppi", "progress_file", "blood_object", "lung_object"))

# Function to run SCORPION for all cell types
run_scorpion <- function(object, cell_types, output_file) {
    scorpion_output <- list()

    foreach(
        cell_type = cell_types,
        .combine = "c",
        .packages = c("Seurat", "SCORPION", "Matrix")
    ) %dopar% {
        cat(paste(Sys.time(), "- Starting:", cell_type, "\n"), file = progress_file, append = TRUE)

        # Subset data by cell type
        cell_subset <- subset(object, idents = cell_type)

        # Fetch gene expression matrix as sparse matrix
        gex_cell <- as(cell_subset@assays$RNA@counts, "dgCMatrix")

        # Create list for SCORPION input
        scorpion_input <- list(gex = gex_cell, tf = tf, ppi = ppi)

        # Run SCORPION
        result <- SCORPION::scorpion(
            tfMotifs = scorpion_input$tf,
            gexMatrix = scorpion_input$gex,
            ppiNet = scorpion_input$ppi,
            alphaValue = 0.1
        )

        if (is.null(result)) {
            cat("SCORPION failed for", cell_type, "\n", file = progress_file, append = TRUE)
            return(NULL)
        }

        save(result, file = sub("output/networks/temp", paste0("output/networks/temp/", cell_type, "_"), output_file))

        # Clear memory
        rm(cell_subset, gex_cell, scorpion_input)
        gc()

        cat(paste(Sys.time(), "- Finished:", cell_type, "\n"), file = progress_file, append = TRUE)

        # Return result
        list(cell_type = result)
    }

    # Save output
    save(scorpion_output, file = output_file)
    cat("Saved SCORPION output -", output_file, "\n", file = progress_file, append = TRUE)
}

# Run SCORPION for blood
run_scorpion(blood_object, c("cd4_positive_t_cell", "cd8_positive_t_cell", "monocyte", "platelet"), "Rdata/bloodScorpionOutput.Rdata")

# Run SCORPION for lung
# run_scorpion(lung_object, c("cd4_positive_t_cell", "cd8_positive_t_cell", "monocyte", "type_ii_pneumocyte"), "Rdata/lungScorpionOutput.Rdata")

# Stop parallel processing
stopCluster(cl)
cat("All tasks completed at", Sys.time(), "\n", file = progress_file, append = TRUE)
