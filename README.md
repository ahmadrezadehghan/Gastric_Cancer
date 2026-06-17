# 🧬 Gastric Cancer Prognostic Gene Signature Pipeline

[![R](https://img.shields.io/badge/R-≥4.0-blue.svg)](https://www.r-project.org/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

**A complete, reproducible pipeline for developing and validating a prognostic gene expression signature for gastric cancer, using multiple public cohorts.**

---

## 📖 Table of Contents

- [Overview](#overview)
- [Features](#features)
- [Pipeline Steps](#pipeline-steps)
  - [1. Preprocessing (no data leakage)](#1-preprocessing-no-data-leakage)
  - [2. Diagnostic Validation](#2-diagnostic-validation)
  - [3. Survival Modelling](#3-survival-modelling)
  - [4. Country & Region-Specific Analysis](#4-country--region-specific-analysis)
  - [5. Network Analysis](#5-network-analysis)
  - [6. Publication-Ready Visualizations](#6-publication-ready-visualizations)
- [Requirements & Installation](#requirements--installation)
- [Input Data Format](#input-data-format)
- [Usage](#usage)
- [Outputs](#outputs)
- [Citation](#citation)
- [License](#license)

---

## Overview

This repository provides a fully documented R pipeline for the discovery and validation of a prognostic gene signature in gastric cancer. It integrates RNA‑seq and microarray data from multiple cohorts, handles batch effects without data leakage, builds an elastic‑net Cox model, and performs comprehensive country‑ and region‑specific survival analyses.

The pipeline is designed to meet the standards of high‑impact journals (e.g., *Nature Medicine*) by ensuring:
- Strict separation of training and test sets at the study level.
- All normalisation parameters learned exclusively from training data.
- Properly weighted decision curve analysis using IPCW weights.
- Extensive QC and diagnostic checks.

---

## ✨ Features

- **No data leakage** – per‑study stratified split; normalisation parameters are locked from training and applied to test.
- **Multi‑platform support** – handles RNA‑seq (DESeq2 VST) and microarray (quantile normalisation) in one workflow.
- **Batch effect correction** – optional ComBat + fsva, with QC plots before/after.
- **Survival modelling** – elastic‑net Cox regression with single or multiple imputation (MICE), time‑dependent AUC with 95% CIs, and DCA.
- **Geographic & regional analysis** – world maps of sample distribution, country‑level hazard ratios, subregion KM curves, and region‑specific prognostic gene discovery.
- **Network visualisation** – STRING PPI network of top genes, coloured by HR direction.
- **Publication‑ready figures** – clinical overview, model performance, KM curves, AUC, and DCA plots.

---

## Pipeline Steps

All scripts are in R and should be executed in the order below.

### 1. Preprocessing (no data leakage)  
**Script:** `step1.1_preprocessing.R`  
- Reads expression matrix, study mapping, and clinical data.  
- Performs stratified train/test split per study (70%/30%).  
- Normalises training data using **DESeq2’s VST** for RNA‑seq (estimating size factors and dispersions on training data) and quantile normalisation for microarray.  
- Applies the same transformation to test RNA‑seq data using the size factors and dispersion trend learned from training (ensuring zero leakage).  
- Global quantile harmonisation (target from training).  
- Optional batch correction (ComBat + fsva).  
- Feature selection by median variance across studies.  
- Saves `train_data.csv`, `test_data.csv`, metadata, and QC plots.

### 2. Diagnostic Validation  
**Script:** `step1.2_diagnostic.R`  
- Loads normalised matrices from Step 1.  
- Generates PCA of combined train/test, coloured by study and set.  
- Computes batch R² distribution and silhouette analysis on training set.  
- Produces density overlap plots and a JSON summary of metrics.

### 3. Survival Modelling  
**Script:** `step2_survival.R`  
- Loads clinical data (OS, stage, age, gender).  
- Selects top variable genes and standardises using training statistics.  
- Builds elastic‑net Cox model (α = 0.9) with **single** or **multiple imputation** (MICE).  
- Evaluates on test set: C‑index, time‑dependent AUC (1–5 years) with bootstrap CIs.  
- Compares with a clinical‑only model.  
- Performs DCA for 3‑year survival using correct IPCW weights from training censoring distribution.  
- Outputs risk scores, coefficients, DCA results, and summary JSON.

### 4. Country & Region-Specific Analysis  
**Script:** `step3_country_region.R`  
- Cleans country names and maps sample distribution.  
- Computes country‑level hazard ratios (countries with ≥5 samples).  
- Groups countries into subregions; runs subregion KM and Cox regression.  
- For each subregion (≥30 samples, ≥10 events), performs univariate Cox on top 2,000 variable genes (training only).  
- Identifies significant genes (FDR < 0.05), creates volcano plots and overlap analyses.  
- Performs Hallmark pathway enrichment (via `msigdbr` or built‑in fallback).  
- Generates heatmap of common prognostic genes and pairwise HR scatter plots.

### 5. Network Analysis  
**Script:** `network_analysis_top_hr_genes.R`  
- Reads univariate Cox results for a chosen region.  
- Selects top 150 genes with HR > 1 and top 150 with HR < 1.  
- Maps genes to STRING and retrieves protein‑protein interactions (confidence ≥ 0.4).  
- Builds an igraph network, computes degree, and visualises with `ggplot2`.  
- Exports node/edge tables, unmapped genes list, and a publication‑ready network plot.

### 6. Publication-Ready Visualizations  
**Script:** `visualizations.R`  
- Creates a comprehensive clinical overview figure (10 panels).  
- Generates model performance figures: KM by risk group, risk score distribution, time‑dependent AUC, and model comparison (if available).  
- Saves all figures as high‑resolution PNG/PDF in a dedicated folder.

---

## 🛠 Requirements & Installation

### R (≥ 4.0)

Install required packages from CRAN and Bioconductor:

```r
# CRAN packages
install.packages(c(
  "readxl", "tidyverse", "matrixStats", "jsonlite", "preprocessCore",
  "ggplot2", "RColorBrewer", "survival", "glmnet", "timeROC", "cluster",
  "ggrepel", "rnaturalearth", "rnaturalearthdata", "sf", "viridis",
  "pheatmap", "UpSetR", "patchwork", "httr", "MASS"
))

# Bioconductor packages
if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install(c("DESeq2", "limma", "sva", "STRINGdb", "msigdbr"))
```

> **Note:** `rnaturalearth` may require additional system libraries (e.g., `libudunits2-dev`, `libgdal-dev` on Linux). For macOS, install via Homebrew or follow the package documentation.

---

## 📂 Input Data Format

Place your input files in the project directory (paths editable in each script). The expected structure is:

| File | Description |
|------|-------------|
| `Gene-Set/Gene_Set.xlsx` | Expression matrix: first column = gene identifiers, subsequent columns = sample IDs. |
| `Gene-Set/study.xlsx` | Study mapping: each row contains a study name followed by sample IDs from that study. |
| `Clinical_Data/demo.xlsx` | Clinical features: first column = feature names (e.g., `Overall_survival(Months)`, `Survival status`, `Stage`, `Age`, `Gender`, `Country`). |
| `Annotation/probe_to_gene.csv` | (Optional) Probe‑to‑gene mapping if expression row names are probe IDs. |
| `Gene-Set/frame.xlsx` | Gene symbol ↔ ID mapping (used in region‑specific analysis). |

---

## 🚀 Usage

1. **Clone the repository** and set your working directory to the project root.
2. **Edit the configuration** at the top of each R script:
   - `BASE_DIR` – path to your project folder.
   - `PLATFORM` – choose `"rnaseq"` or `"microarray"` (must match your data).
   - Adjust other parameters (e.g., `TRAIN_FRACTION`, `N_TOP_GENES`, imputation method) as needed.
3. **Run the scripts in order** using R or `Rscript`:

```bash
Rscript step1.1_preprocessing.R
Rscript step1.2_diagnostic.R
Rscript step2_survival.R
Rscript step3_country_region.R
Rscript network_analysis_top_hr_genes.R   # edit REGION variable inside
Rscript visualizations.R
```

> ⚠️ **Important:** The internal train/test split is **not** a substitute for external validation. For high‑impact publications, we strongly recommend validating the final signature on an independent external cohort.

---

## 📤 Outputs

All outputs are organised in subdirectories under `BASE_DIR`:

| Directory | Contents |
|-----------|----------|
| `Prepared_Data/` | Normalised train/test expression matrices, sample metadata, QC plots, diagnostic validation (metrics, figures). |
| `Survival_ElasticNet_DCA/` | Test risk scores, model coefficients, DCA results, summary JSON, and log file. |
| `Country_Region_Analysis/` | Country maps, subregion survival results, region‑specific Cox results, volcano plots, enrichment tables, all‑p‑values summary. |
| `Network_Analysis/` | Network node/edge tables, unmapped genes list, network plots (PNG/PDF), and summary JSON. |
| `Manuscript_Figures/` | Publication‑ready figures: clinical overview, survival performance (KM, risk distribution, AUC, model comparison). |

---

## 📝 Citation

If you use this pipeline in your research, please cite this repository and, where applicable, the original public datasets and R packages used.

---

## 📄 License

This project is released under the [MIT License](LICENSE).

---

## 🙏 Acknowledgements

We thank the developers of all R packages used and the public cohorts that made this analysis possible.

---

**Happy analysing!**  
For questions or issues, please open an issue on GitHub.
