# GRN construction and analysis pipeline
An R pipeline developed for constructing and analyzing gene regulatory networks (GRNs) using scRNA-seq data and SCORPION.\
The Tabula Sapiens v2 dataset was used for this initial construction and analysis of networks.

## 1. Repository structure
```
📁 Pipeline/
    📁 data/
        📁 config/
            📄 config.yaml
        📁 priors/
            📄 ppi_prior_2024.tsv                          # Not tracked due to file size
            📄 motif_prior_names_2024.tsv                  # Not tracked due to file size
        📁 rds/
    📁 output/
        📁 analysis/
            📁 correlation/
            📁 dimreduction/
                📁 umap/
            📁 gsea_dotplots/
                📁 gex/
                📁 indegree/
            📁 linear_regression/
                📁 gex/
                📁 indegree/
        📁 networks/
            📁 final/
            📁 temp/
        📁 preprocessing/
            📁 all/
            📁 plots/
    📁 scripts/
        📄 00_download_TS_v2.sh
        📄 001_installPackages.R
        📄 01_preprocessing.R
        📄 02_runScorpion.R
        📄 03_analyzeNetworks.R
📄 README.md  
📄 .gitignore
```

## 2. Installation & Dependencies
### 1. Install R (version 4.2.1)
### 2. Install required packages by running
```bash
Rscript Pipeline/scripts/001_installPackages.R
```
### 3. Download input files
The scRNA-seq data used for this specific analysis can be downloaded by running
```
bash Pipeline/scripts/00_download_TS_v2.sh
```
The .rds files downloaded will be saved in `Pipeline/data/rds/`

The Seurat objects acquired from Tabula Sapiens contain mistakes and unwanted celltype annotations and needed correction, that's why a config file was used for this analysis ([Pipeline/data/config/config.yaml](https://github.com/siepvan020/reference-networks/blob/1165731dea81a9494c16b94c7939ad3f82167581/Pipeline/data/config/config.yaml)). Also celltypes that do not need correction are defined here. The YAML file is structured as following:
```
{tissue}:                        # Tissue name used in R scripts - lowercase and no whitespaces!!!
  file: TS_v2_{Tissue}.rds       # Input filename (will be automatically correct after running Pipeline/scripts/00_download_TS_v2.sh
  cells:
    {corrected_celltype}:        # Final orrected celltype annotation - lowercase and no whitespaces!!!
      - Cell Type 1A             # Incorrect/redundant celltype annotation to be corrected
      - Cell type 1B             # Incorrect/redundant celltype annotation to be corrected
    {celltype}:                  # Celltype that does not need annotation correction
```


The prior support networks are not included in the repository, these can be accessed via https://zenodo.org/records/15063580 \
The prior networks used for this analysis can be accessed via:
- https://zenodo.org/records/15063580/files/motif_prior_names_2024.tsv
- https://zenodo.org/records/15063580/files/ppi_prior_2024.tsv

The pipeline assumes these prior networks are stored in `Pipeline/data/priors/`

### 4. Running the Pipeline
Run scripts in order:
```r
source("Pipeline/scripts/01_preprocessing.R")
source("Pipeline/scripts/02_runscorpion.R")
source("Pipeline/scripts/03_analyzenetworks.R")
```
### 5. Output files
#### The script **01_preprocessing.R** produces the following files:
- Pipeline/output/preprocessing/all/{tissue}_prepped.rds
- Pipeline/output/preprocessing/plots/all_tissues_batch_effect_umap.pdf

Where all the configured tissues will have one preprocessed (prepped) RDS file.

#### The script **02_runScorpion.R** produces the following files:
- Pipeline/output/networks/final/{tissue}ScorpionOutput.Rdata
- Pipeline/output/networks/temp/{celltype}_{tissue}ScorpionOutput.Rdata

Where every configurated tissue's RData object in `final/` contains a list of celltype-specific networks.
The RData files in the `temp/` folder contain one network per file, and are stored to save the network if anything unexpected happens when saving the final object.

#### The script **03_analyzeNetworks** produces files in the following directories (too many to list all files):
- Pipeline/output/analysis/dimreduction/umap/
- Pipeline/output/analysis/correlation/
- Pipeline/output/analysis/linear_regression/gex
- Pipeline/output/analysis/linear_regression/indegree/
- Pipeline/output/analysis/gsea_dotplots/gex
- Pipeline/output/analysis/gsea_dotplots/indegree
