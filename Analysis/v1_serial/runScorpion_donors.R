#!/usr/bin/env Rscript

# Define log file
log_file <- "/div/pythagoras/u1/siepv/siep/Analysis/logfile_donors.txt"

# # Open the log file connection
log_con <- file(log_file, open = "wt")

# Start logging
sink(file = log_con, append = TRUE) # Redirect standard output
sink(file = log_con, append = TRUE, type = "message") # Redirect error messages

cat("#######################################\nScript started at", Sys.time(), "\n")


# Set CRAN mirror
options(repos = c(CRAN = "https://cloud.r-project.org"))

# Load required packages
if (!requireNamespace("remotes", quietly = TRUE)) install.packages("remotes")
if (!requireNamespace("SeuratDisk", quietly = TRUE)) remotes::install_github("mojaveazure/seurat-disk")
if (!requireNamespace("SCORPION", quietly = TRUE)) install.packages("SCORPION")

# Load required packages
library(SeuratDisk)
library(SCORPION)
library(Seurat)

# Set working directory
setwd("/div/pythagoras/u1/siepv/siep/Analysis")

# Load Seurat objects
readRDS("Rdata/blood_object.rds") -> blood_object
readRDS("Rdata/lung_object.rds") -> lung_object

# Read priors
ppi <- read.delim("data/priors/ppi_prior_2024.tsv", header = FALSE, sep = "\t")
tf <- read.delim("data/priors/motif_prior_names_2024.tsv", header = FALSE, sep = "\t")

# Run SCORPION for blood
bloodScorpionOutput <- list()
blood_cell_types <- list("cd4_positive_t_cell", "cd8_positive_t_cell", "monocyte", "platelet")
blood_donors <- unique(blood_object$donor)

for (cell_type in blood_cell_types) {
    for (donor in blood_donors) {
        cat("--------------------", cell_type, " - Donor:", donor, "--------------------\n")

        # Subset data by cell type and donor
        blood_cell_subset <- subset(blood_object, idents = cell_type, subset = donor == donor)

        # Set cutoff for > 100 cells
        if (ncol(blood_cell_subset) < 100) {
            cat("Skipped:", cell_type, "- Donor:", donor, "- Not enough cells (", ncol(blood_cell_subset), ")\n")
            next # Skip this donor-cell type combination
        }

        # Fetch gene expression matrix
        gex_cell <- blood_cell_subset@assays$decontXcounts@counts

        # Create list for SCORPION input
        bloodScorpionInput <- list(gex = gex_cell, tf = tf, ppi = ppi)

        # Run SCORPION and store the result
        bloodScorpionOutput[[paste(cell_type, donor, sep = "_")]] <- SCORPION::scorpion(
            tfMotifs = bloodScorpionInput$tf,
            gexMatrix = bloodScorpionInput$gex,
            ppiNet = bloodScorpionInput$ppi,
            alphaValue = 0.1
        )
    }
}

save(bloodScorpionOutput, "Rdata/bloodScorpionOutput_donors.Rdata")
cat("Saved SCORPION output - Blood tissue")

# Run SCORPION for lung
lungScorpionOutput <- list()
lung_cell_types <- list("cd4_positive_t_cell", "cd8_positive_t_cell", "monocyte", "type_ii_pneumocyte")
lung_donors <- unique(lung_object$donor)

for (cell_type in lung_cell_types) {
    for (donor in lung_donors) {
        cat("--------------------", cell_type, " - Donor:", donor, "--------------------\n")

        # Subset data by cell type and donor
        lung_cell_subset <- subset(lung_object, idents = cell_type, subset = donor == donor)

        # Set cutoff for > 100 cells
        if (ncol(lung_cell_subset) < 100) {
            cat("Skipped:", cell_type, "- Donor:", donor, "- Not enough cells (", ncol(lung_cell_subset), ")\n")
            next # Skip this donor-cell type combination
        }

        # Fetch gene expression matrix
        gex_cell <- lung_cell_subset@assays$decontXcounts@counts

        # Create list for SCORPION input
        lungScorpionInput <- list(gex = gex_cell, tf = tf, ppi = ppi)

        # Run SCORPION and store the result
        lungScorpionOutput[[paste(cell_type, donor, sep = "_")]] <- SCORPION::scorpion(
            tfMotifs = lungScorpionInput$tf,
            gexMatrix = lungScorpionInput$gex,
            ppiNet = lungScorpionInput$ppi,
            alphaValue = 0.1
        )
    }
}

save(lungScorpionOutput, "Rdata/lungScorpionOutput_donors.Rdata")
cat("Saved SCORPION output - Lung tissue")

# Reset sinks and close the log file connection
while (sink.number() > 0) sink()

sink(type = "message")
