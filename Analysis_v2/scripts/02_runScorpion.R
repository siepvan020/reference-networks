#!/usr/bin/env Rscript

#### 1. Setup ####

progress_file <- "/div/pythagoras/u1/siepv/siep/Analysis_v2/output/log/network_progress.log"

# Start fresh log or append
cat("Starting script at", format(Sys.time()), "\n", file = progress_file)

# Set CRAN mirror
options(repos = c(CRAN = "https://cran.r-project.org"))

# Install required packages
# if (!requireNamespace("remotes", quietly = TRUE)) install.packages("remotes")
# if (!requireNamespace("Matrix", quietly = TRUE)) remotes::install_version("Matrix", version = "1.6.4")
# install.packages("SeuratObject")
# if (!requireNamespace("SeuratObject", quietly = TRUE)) remotes::install_version("SeuratObject", "5.2.0", repos = c("https://satijalab.r-universe.dev", getOption("repos")))
# if (!requireNamespace("Seurat", quietly = TRUE)) remotes::install_version("Seurat", "5.2.0", repos = c("https://satijalab.r-universe.dev", getOption("repos")))
if (!requireNamespace("Seurat", quietly = TRUE)) install.packages("Seurat")
if (!requireNamespace("SCORPION", quietly = TRUE)) install.packages("SCORPION")
if (!requireNamespace("doParallel", quietly = TRUE)) install.packages("doParallel")
if (!requireNamespace("Matrix", quietly = TRUE)) install.packages("Matrix")

# Load required packages
library(Seurat)
library(SeuratObject)
library(SCORPION)
library(doParallel)
library(Matrix)


# Set working directory
setwd("/div/pythagoras/u1/siepv/siep/Analysis_v2")


#### 2. Load data ####

# Function to load Seurat objects and set cell type identities
load_seurat_object <- function(file_path) {
    object <- readRDS(file_path)
    Idents(object) <- "cell_type"
    return(object)
}

# Load Seurat objects
blood_object <- load_seurat_object("output/preprocessing/blood_prepped.rds")
lung_object <- load_seurat_object("output/preprocessing/lung_prepped.rds")
# fat_object <- load_seurat_object("output/preprocessing/fat_filter.rds")
# kidney_object <- load_seurat_object("output/preprocessing/kidney_filter.rds")
# liver_object <- load_seurat_object("output/preprocessing/liver_filter.rds")

cat("- Loaded Seurat objects -", format(Sys.time()), "\n", file = progress_file, append = TRUE)

# Read priors
ppi <- read.delim("data/priors/ppi_prior_2024.tsv", header = FALSE, sep = "\t")
tf <- read.delim("data/priors/motif_prior_names_2024.tsv", header = FALSE, sep = "\t")
cat("- Loaded priors -", format(Sys.time()), "\n", file = progress_file, append = TRUE)

# Set up parallel processing
num_cores <- 4
cl <- makeCluster(num_cores)
registerDoParallel(cl)

#### 3. Run SCORPION for each celltype####
run_scorpion <- function(tissue, object, cell_types, output_file) {
    scorpion_output <- list()

    scorpion_output <- foreach(
        cell_type = cell_types,
        .combine = "c",
        .packages = c("Seurat", "SCORPION", "Matrix"),
        .errorhandling = "pass",
        .verbose = TRUE,
        .export = c("tf", "ppi", "progress_file", "object")
    ) %dopar% {
        cat("\n- Starting:", cell_type, " -", format(Sys.time()), "\n", file = progress_file, append = TRUE)

        # Subset data by cell type
        cell_subset <- subset(object, idents = cell_type)

        # Fetch gene expression matrix as sparse matrix
        gex_cell <- as(cell_subset@assays$RNA@counts, "dgCMatrix")

        # Create list for SCORPION input
        scorpion_input <- list(gex = gex_cell, tf = tf, ppi = ppi)

        #### - Run SCORPION algorithm ####
        result <- tryCatch(
            {
                SCORPION::scorpion(
                    tfMotifs = scorpion_input$tf,
                    gexMatrix = scorpion_input$gex,
                    ppiNet = scorpion_input$ppi,
                    alphaValue = 0.1 # ,
                    # outNet = c("regulatory")
                )
            },
            error = function(e) {
                cat("\nSCORPION failed for", cell_type, format(Sys.time()), "with error:", e$message, "\n", file = progress_file, append = TRUE)
                return(NULL)
            }
        )

        if (is.null(result)) {
            return(NULL)
        }
        if ("regNet" %in% names(result)) {
            result <- result$regNet
        }

        #### - Save and return result ####
        temp_file <- file.path("output/networks/temp/", paste0(cell_type, "_", basename(output_file)))
        save(result, file = temp_file)

        # Clear memory
        rm(cell_subset, gex_cell, scorpion_input)
        gc()

        cat("\n- Finished:", cell_type, " -", format(Sys.time()), "\n", file = progress_file, append = TRUE)

        # Return result
        setNames(list(result), cell_type)
    }

    # Dynamically assign the scorpion_output variable name based on the tissue parameter
    assign(paste0(tissue, "_output"), scorpion_output, envir = .GlobalEnv)

    # Save the dynamically named variable to the output file
    save(list = paste0(tissue, "_output"), file = output_file)
    cat("Saved SCORPION output -", output_file, "\n", file = progress_file, append = TRUE)
}

# Run SCORPION for blood
run_scorpion("blood", blood_object, c("cd4_positive_t_cell", "cd8_positive_t_cell", "monocyte", "platelet"), "output/networks/final/bloodScorpionOutput4.Rdata")

# Run SCORPION for lung
run_scorpion("lung", lung_object, c("cd4_positive_t_cell", "cd8_positive_t_cell", "monocyte", "type_ii_pneumocyte"), "output/networks/final/lungScorpionOutput4.Rdata")

# # Run SCORPION for fat
# run_scorpion("fat", fat_object, c("cd4_positive_t_cell", "cd8_positive_t_cell", "monocyte", "adipose_stem_cell"), "output/networks/final/fatScorpionOutput.Rdata")

# # Run SCORPION for kidney
# run_scorpion("kidney", kidney_object, c("cd4_positive_t_cell", "cd8_positive_t_cell", "monocyte", "kidney_epithelial_cell"), "output/networks/final/kidneyScorpionOutput.Rdata")

# # Run SCORPION for liver
# run_scorpion("liver", liver_object, c("cd4_positive_t_cell", "cd8_positive_t_cell", "monocyte", "hepatocyte"), "output/networks/final/liverScorpionOutput.Rdata")

# Stop parallel processing
stopCluster(cl)
cat("All tasks completed at", format(Sys.time()), "\n", file = progress_file, append = TRUE)
