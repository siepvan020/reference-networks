#!/usr/bin/env Rscript

#### 1. Setup ####
# Helper function
`%||%` <- function(a, b) if (!is.null(a)) a else b

# Ensure that remotes is installed
if (!requireNamespace("remotes", quietly = TRUE)) {
  install.packages("remotes")
}


#### 2. Installment functions ####

install_package <- function(pkg, version = NULL, repo = "https://cloud.r-project.org") {
  installed <- requireNamespace(pkg, quietly = TRUE)
  current_version <- if (installed) as.character(utils::packageVersion(pkg)) else NA

  if (is.null(version) && !installed) {
    install.packages(pkg, repos = repo, lib = lib_path)
  }

  if (!is.null(version) && (!installed || current_version != version)) {
    remotes::install_version(pkg, version = version, repos = repo, lib = lib_path)
  }
}

install_from_url <- function(pkg, version, url) {
  installed <- requireNamespace(pkg, quietly = TRUE)
  current_version <- if (installed) as.character(utils::packageVersion(pkg)) else NA

  if (!installed || current_version != version) {
    remotes::install_url(url, lib = lib_path)
  }
}


#### 3. Preinstall Seurat dependencies ####

# Install Seurat dependencies if needed
install_from_url("spatstat.explore", "3.3-4", "https://cran.r-project.org/src/contrib/Archive/spatstat.explore/spatstat.explore_3.3-4.tar.gz")
install_from_url("spatstat.geom", "3.3-4", "https://cran.r-project.org/src/contrib/Archive/spatstat.geom/spatstat.geom_3.3-4.tar.gz")
# install_url_if_needed("sp", "2.1-4", "https://cran.r-project.org/src/contrib/Archive/sp/sp_2.1-4.tar.gz")

install_package("SeuratObject", "4.1.4", repo = "https://satijalab.r-universe.dev")


#### 4. Define remaining packages with their versions ####

# List of packages (names = package names, values = version or NULL for latest)
packages <- list(
    "BiocParallel" = NULL,
    "biomaRt" = "2.52.0",
    "doParallel" = "1.0.17",
    "dplyr" = "1.1.4",
    "fgsea" = "1.22.0",
    "genefilter" = NULL,
    "ggplot2" = "3.5.1",
    "ggrepel" = "0.9.6",
    "glue" = NULL,
    "iterators" = "1.0.14",
    "Matrix" = "1.6-4",
    "msigdbr" = "10.0.1",
    "patchwork" = "1.3.0",
    "purrr" = "1.0.4",
    "sp" = "2.2-0",
    "SCORPION" = "1.0.2",
    "stringr" = "1.5.1",
    "sva" = "3.44.0",
    "tidyr" = "1.3.1",
    "uwot" = NULL,
    "yaml" = "2.3.10"
)

# Custom repos if applicable
custom_repos <- list(
  "sp" = "https://edzer.r-universe.dev"
)


#### 5. Loop through packages and install them ####

for (pkg in names(packages)) {
  repo <- custom_repos[[pkg]] %||% "https://cloud.r-project.org"
  install_package(pkg, packages[[pkg]], repo)
}


#### 6. Install Seurat seperately without dependencies ####
if (!requireNamespace("Seurat", quietly = TRUE) ||
    as.character(utils::packageVersion("Seurat")) != "4.4.0") {
  remotes::install_version("Seurat", version = "4.4.0", repos = "https://satijalab.r-universe.dev", dependencies = TRUE, lib = "/tmp/test_lib")
}