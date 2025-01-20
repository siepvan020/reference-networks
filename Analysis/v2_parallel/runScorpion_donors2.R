#!/usr/bin/env Rscript

# Define log file
# log_file <- "/div/pythagoras/u1/siepv/siep/Analysis/logfile_donors_parallel.txt"

# # Open the log file connection
# log_con <- file(log_file, open = "wt")

# # Start logging
# sink(file = log_con, append = TRUE) # Redirect standard output
# sink(file = log_con, append = TRUE, type = "message") # Redirect error messages

# cat("#######################################\nScript started at", Sys.time(), "\n")

progress_file <- "/div/pythagoras/u1/siepv/siep/Analysis/v2_parallel/progress_donors.log"

# Start fresh log or append
cat("Starting script at", Sys.time(), "\n", file = progress_file)

# Set CRAN mirror
options(repos = c(CRAN = "https://cloud.r-project.org"))

# Load required packages
if (!requireNamespace("remotes", quietly = TRUE)) install.packages("remotes")
if (!requireNamespace("SeuratDisk", quietly = TRUE)) remotes::install_github("mojaveazure/seurat-disk")
if (!requireNamespace("SCORPION", quietly = TRUE)) install.packages("SCORPION")
if (!requireNamespace("doParallel", quietly = TRUE)) install.packages("doParallel")

# Load required packages
library(SeuratDisk)
library(SCORPION)
library(Seurat)
library(doParallel)

# Set working directory
setwd("/div/pythagoras/u1/siepv/siep/Analysis")

# Load Seurat objects
blood_object <- readRDS("Rdata/blood_object.rds")
lung_object <- readRDS("Rdata/lung_object.rds")
cat("- Loaded Seurat objects - \n", file = progress_file, append = TRUE)

# Read priors
ppi <- read.delim("data/priors/ppi_prior_2024.tsv", header = FALSE, sep = "\t")
tf <- read.delim("data/priors/motif_prior_names_2024.tsv", header = FALSE, sep = "\t")
cat("- Loaded priors - \n", file = progress_file, append = TRUE)

# Set up parallel processing
num_cores <- 4
cl <- makeCluster(num_cores)
registerDoParallel(cl)
clusterExport(cl, varlist = c("tf", "ppi", "progress_file", "blood_object", "lung_object"))

# Function to run SCORPION for all donors and cell types
run_scorpion_donor <- function(object, cell_types, donors, output_file) {
    scorpion_output <- list()

    foreach(
        cell_type = cell_types,
        .combine = "c",
        .packages = c("Seurat", "SCORPION"),
    ) %:% foreach(
        donor = donors,
        .combine = "c"
    ) %dopar% {
        cat(paste(Sys.time(), "- Starting:", cell_type, "- Donor:", donor, "\n"), file = progress_file, append = TRUE)

        # Subset data by cell type and donor
        cell_subset <- subset(object, idents = cell_type, subset = donor == donor)

        # Check if subset meets cell cutoff
        if (ncol(cell_subset) < 100) {
            cat(paste("Skipped:", cell_type, "- Donor:", donor, "- Not enough cells (", ncol(cell_subset), ")\n"), file = progress_file, append = TRUE)
            return(NULL) # Skip this combination
        }

        # Fetch gene expression matrix
        gex_cell <- cell_subset@assays$decontXcounts@counts

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

        save(result, file = sub("Rdata/", paste0("Rdata/temp_", cell_type, "_"), output_file))

        # Clear memory
        rm(cell_subset, gex_cell, scorpion_input)
        gc()

        cat(paste(Sys.time(), "- Finished:", cell_type, "- Donor:", donor, "\n"), file = progress_file, append = TRUE)

        # Return result
        setNames(list(result), paste(cell_type, donor, sep = "_"))
    }

    # Save output
    save(scorpion_output, file = output_file)
    cat(paste(Sys.time(), "- Saved SCORPION output to:", output_file, "\n"), file = progress_file, append = TRUE)
}

# Run SCORPION for blood
run_scorpion_donor(blood_object, c("cd4_positive_t_cell", "cd8_positive_t_cell", "monocyte", "platelet"), unique(blood_object$donor), "Rdata/bloodScorpionOutput_donors.Rdata")

# Run SCORPION for lung
run_scorpion_donor(lung_object, c("cd4_positive_t_cell", "cd8_positive_t_cell", "monocyte", "type_ii_pneumocyte"), unique(lung_object$donor), "Rdata/lungScorpionOutput_donors.Rdata")

# Stop parallel processing
stopCluster(cl)
cat("All tasks completed at", Sys.time(), "\n", file = progress_file, append = TRUE)

# # Reset sinks and close the log file connection
# while (sink.number() > 0) sink()
# sink(type = "message")
