#!/usr/bin/env Rscript

#### 1. Setup ####

progress_file <- "/div/pythagoras/u1/siepv/siep/Pipeline/output/log/networks.log"

# Start fresh log or append
cat("Starting script at", format(Sys.time()), "\n", file = progress_file)

# Load required packages
library(Seurat)
library(SeuratObject)
library(SCORPION)
library(doParallel)
library(Matrix)

library(yaml)
library(glue)
library(purrr)

# Set working directory
setwd("/div/pythagoras/u1/siepv/siep/Pipeline")

# Load config file
config <- yaml::read_yaml("data/config/config.yaml")


#### 2. Define helper functions ####

# Function to load Seurat objects and set cell type identities
load_seurat_object <- function(file_path) {
    object <- readRDS(file_path)
    Idents(object) <- "cell_type"
    return(object)
}

# Function to run SCORPION for a given tissue and cell types
run_scorpion <- function(object, cell_type, tissue, outfile) {

  cell_subset <- subset(object, idents = cell_type)
  gex         <- as(cell_subset@assays$RNA@counts, "dgCMatrix")

  res <- tryCatch(
    SCORPION::scorpion(
      tfMotifs  = tf,
      gexMatrix = gex,
      ppiNet    = ppi,
      alphaValue = 0.1
    )$regNet,
    error = function(e) {
      cat("SCORPION failed for", tissue, cell_type, ":", e$message, "\n",
          file = progress_file, append = TRUE)
      NULL
    }
  )
  
  temp_path <- file.path("output/networks/temp/", paste0(cell_type, "_", basename(outfile)))
  save(res, file = temp_path)

  rm(cell_subset, gex)
  gc()
  res
}


#### 3. Read priors and set up parallel processing ####

# Read priors
ppi <- read.delim("data/priors/ppi_prior_2024.tsv", header = FALSE, sep = "\t")
tf <- read.delim("data/priors/motif_prior_names_2024.tsv", header = FALSE, sep = "\t")
cat("- Loaded priors -", format(Sys.time()), "\n", file = progress_file, append = TRUE)

# Set up parallel processing
num_cores <- 2
cl <- makeCluster(num_cores)
registerDoParallel(cl)


#### 4. Run SCORPION for each tissue in the config file ####

process_tissue <- function(spec, tissue) {

  cat("\n=== Tissue:", tissue, " ===", format(Sys.time()), "\n",
      file = progress_file, append = TRUE)

  # 4.1  load the pre-processed Seurat object
  obj_path <- glue("output/preprocessing/all/{tissue}_prepped.rds")
  object   <- load_seurat_object(obj_path)

  cell_types <- names(spec$cells)
  outfile    <- glue("output/networks/final/{tissue}ScorpionOutput.Rdata")

  # 4.2  parallel loop over cell-types
  results <- foreach(
      ct        = cell_types,
      .combine  = function(a, b) c(a, b),     # keep names
      .packages = c("Seurat", "SCORPION", "Matrix"),
      .export   = c("tf", "ppi", "progress_file", "run_scorpion", "object", "outfile", "tissue")
    ) %dopar% {

      cat("- Starting:", ct, format(Sys.time()), "\n",
          file = progress_file, append = TRUE)

      res <- run_scorpion(object, ct, tissue, outfile)

      cat("- Finished:", ct, format(Sys.time()), "\n",
          file = progress_file, append = TRUE)

      setNames(list(res), ct)
    }

  # 4.3  save once per tissue
  assign(paste0(tissue, "_output"), results, envir = .GlobalEnv)
  save(list = paste0(tissue, "_output"), file = outfile)

  cat("Saved:", outfile, "\n", file = progress_file, append = TRUE)

  rm(object, results)
  gc()
}

purrr::iwalk(config, process_tissue)

# Stop parallel processing
stopCluster(cl)
cat("All tasks completed at", format(Sys.time()), "\n", file = progress_file, append = TRUE)
