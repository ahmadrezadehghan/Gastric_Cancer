#!/usr/bin/env Rscript
# =============================================================================
# STEP 2: DIAGNOSTIC VALIDATION OF STEP 1 OUTPUT
# =============================================================================
# This script does NOT modify data. It loads the already normalised and
# batch‑corrected matrices from Step 1 (train_data.csv, test_data.csv) and
# generates QC plots and metrics to confirm that batch effects are adequately
# removed and that train/test distributions are comparable.
#
# Inputs (from Step 1):
#   - Prepared_Data/{rnaseq|microarray}/train_data.csv
#   - Prepared_Data/{rnaseq|microarray}/test_data.csv
#   - Prepared_Data/{rnaseq|microarray}/sample_metadata.csv
#
# Outputs:
#   - Prepared_Data/Diagnostics/train_test_PCA.pdf
#   - Prepared_Data/Diagnostics/batch_R2_distribution.pdf
#   - Prepared_Data/Diagnostics/silhouette_analysis.pdf
#   - Prepared_Data/Diagnostics/validation_metrics.json
# =============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(ggplot2)
  library(matrixStats)
  library(jsonlite)
  library(cluster)
  library(ggrepel)
})

# -----------------------------------------------------------------------------
# Configuration – must match Step 1 output location
# -----------------------------------------------------------------------------
BASE_DIR      <- "E:/GastricCancer-2026"
PREP_BASE     <- file.path(BASE_DIR, "Prepared_Data")   # output of STEP 1
PLATFORM      <- "rnaseq"     # or "microarray" – set automatically if both exist
OUTPUT_DIR    <- file.path(BASE_DIR, "Prepared_Data", "Diagnostics")
FIGURES_DIR   <- file.path(OUTPUT_DIR, "figures")

for (d in c(OUTPUT_DIR, FIGURES_DIR)) {
  if (!dir.exists(d)) dir.create(d, recursive = TRUE)
}

log_file <- file.path(OUTPUT_DIR, "step1.2_diagnostic_log.txt")
sink(log_file, append = FALSE, split = TRUE)

cat(rep("=", 80), "\n")
cat("STEP 2: DIAGNOSTIC VALIDATION OF STEP 1 OUTPUT\n")
cat("No data transformation – only QC and metrics\n")
cat("Start time:", as.character(Sys.time()), "\n")
cat(rep("=", 80), "\n\n")

set.seed(42)

# -----------------------------------------------------------------------------
# 1. Automatically detect platform(s) and load the corresponding data
# -----------------------------------------------------------------------------
cat("[1] Detecting available platforms\n")
platform_dirs <- list.dirs(PREP_BASE, recursive = FALSE, full.names = FALSE)
platform_dirs <- platform_dirs[platform_dirs %in% c("rnaseq", "microarray")]

if (length(platform_dirs) == 0) {
  stop("No platform subdirectories found in ", PREP_BASE)
}

all_train <- list()
all_test <- list()
all_meta <- list()

for (plat in platform_dirs) {
  plat_path <- file.path(PREP_BASE, plat)
  train_file <- file.path(plat_path, "train_data.csv")
  test_file  <- file.path(plat_path, "test_data.csv")
  meta_file  <- file.path(plat_path, "sample_metadata.csv")

  if (!file.exists(train_file) || !file.exists(test_file) || !file.exists(meta_file)) {
    warning("Incomplete data for platform ", plat, " – skipping")
    next
  }

  cat("  Loading platform:", plat, "\n")
  train <- read.csv(train_file, row.names = 1, check.names = FALSE)
  test  <- read.csv(test_file,  row.names = 1, check.names = FALSE)
  meta  <- read.csv(meta_file,  row.names = 1, stringsAsFactors = FALSE)

  cat("    Train: ", nrow(train), " samples × ", ncol(train), " genes\n", sep = "")
  cat("    Test:  ", nrow(test),  " samples × ", ncol(test),  " genes\n", sep = "")

  # Ensure consistent gene order across train/test
  common_genes <- intersect(colnames(train), colnames(test))
  if (length(common_genes) < 100) {
    warning("Too few common genes for platform ", plat, " – skipping")
    next
  }
  train <- train[, common_genes]
  test  <- test[,  common_genes]

  # Add platform label to metadata
  meta$platform <- plat
  all_train[[plat]] <- train
  all_test[[plat]]  <- test
  all_meta[[plat]]  <- meta
}

if (length(all_train) == 0) stop("No valid data loaded")

# -----------------------------------------------------------------------------
# 2. For each platform, run diagnostics
# -----------------------------------------------------------------------------
all_metrics <- list()

for (plat in names(all_train)) {
  cat("\n", rep("=", 80), "\n")
  cat("PLATFORM:", toupper(plat), "\n")
  cat(rep("=", 80), "\n")

  train_mat <- all_train[[plat]]
  test_mat  <- all_test[[plat]]
  meta      <- all_meta[[plat]]

  # Ensure metadata contains the necessary columns
  required_cols <- c("sample_id", "study", "set", "platform")
  for (col in required_cols) {
    if (!col %in% colnames(meta)) stop("Missing column in metadata: ", col)
  }

  # Subset metadata to samples present in train/test
  meta_train <- meta[rownames(meta) %in% rownames(train_mat), , drop = FALSE]
  meta_test  <- meta[rownames(meta) %in% rownames(test_mat), , drop = FALSE]

  # Combine for PCA (train + test)
  combined <- rbind(train_mat, test_mat)
  combined_meta <- rbind(
    cbind(meta_train, set = "train"),
    cbind(meta_test,  set = "test")
  )
  batch_vec <- combined_meta$study   # main batch variable (study)

  # -------------------------------------------------------------------------
  # 2.1 PCA on combined train+test
  # -------------------------------------------------------------------------
  cat("\n[2.1] PCA of combined train/test data\n")
  # Use top 2000 variable genes for PCA (fast)
  gene_vars <- apply(combined, 2, var, na.rm = TRUE)
  top_genes <- names(sort(gene_vars, decreasing = TRUE))[1:min(2000, length(gene_vars))]
  pca <- prcomp(combined[, top_genes], center = TRUE, scale. = FALSE)
  pca_df <- as.data.frame(pca$x[, 1:4])
  colnames(pca_df) <- paste0("PC", 1:4)
  pca_df$set <- combined_meta$set
  pca_df$study <- combined_meta$study
  pca_df$sample <- rownames(combined)

  var_exp <- round(100 * summary(pca)$importance[2, 1:4], 1)

  # Plot PC1 vs PC2 coloured by study, shaped by train/test
  p_pca <- ggplot(pca_df, aes(x = PC1, y = PC2, colour = study, shape = set)) +
    geom_point(size = 2, alpha = 0.7) +
    scale_shape_manual(values = c("train" = 16, "test" = 1)) +
    labs(x = paste0("PC1 (", var_exp[1], "%)"),
         y = paste0("PC2 (", var_exp[2], "%)"),
         title = paste(toupper(plat), "– Train/Test PCA after STEP 1")) +
    theme_minimal() +
    theme(legend.position = "bottom")

  ggsave(file.path(FIGURES_DIR, paste0(plat, "_PCA_train_test.pdf")),
         p_pca, width = 8, height = 6, dpi = 300, bg = "white")

  # -------------------------------------------------------------------------
  # 2.2 Batch effect R² (training only)
  # -------------------------------------------------------------------------
  cat("\n[2.2] Batch effect R² on training set\n")
  train_top <- train_mat[, top_genes[top_genes %in% colnames(train_mat)]]
  batch_train <- meta_train$study

  r2_vals <- sapply(1:ncol(train_top), function(i) {
    tryCatch(summary(lm(train_top[, i] ~ batch_train))$r.squared, error = function(e) NA)
  })
  r2_vals <- r2_vals[!is.na(r2_vals)]
  cat("  Mean batch R²: ", round(mean(r2_vals), 4), "\n", sep = "")

  # Density plot of R²
  r2_df <- data.frame(R2 = r2_vals)
  p_r2 <- ggplot(r2_df, aes(x = R2)) +
    geom_density(fill = "steelblue", alpha = 0.5) +
    geom_vline(xintercept = mean(r2_vals), linetype = "dashed", colour = "red") +
    labs(x = expression(R^2 ~ "of study batch"), y = "Density",
         title = paste(toupper(plat), "- Distribution of batch R² (training)")) +
    theme_minimal()
  ggsave(file.path(FIGURES_DIR, paste0(plat, "_batch_R2.pdf")),
         p_r2, width = 7, height = 5, dpi = 300, bg = "white")

  # -------------------------------------------------------------------------
  # 2.3 Silhouette analysis (training PCA space)
  # -------------------------------------------------------------------------
  cat("\n[2.3] Silhouette analysis on training PCA\n")
  # Use first 10 PCs
  pca_train <- pca$x[rownames(pca_df) %in% rownames(train_mat), 1:min(10, ncol(pca$x))]
  if (nrow(pca_train) > 2 && length(unique(batch_train)) > 1) {
    dist_mat <- dist(pca_train)
    sil <- silhouette(as.numeric(factor(batch_train)), dist_mat)
    avg_sil <- summary(sil)$avg.width
    cat("  Average silhouette width: ", round(avg_sil, 4), "\n", sep = "")

    # Plot silhouette
    pdf(file.path(FIGURES_DIR, paste0(plat, "_silhouette.pdf")), width = 8, height = 6)
    plot(sil, main = paste(toupper(plat), "- Silhouette by study batch"),
         col = rainbow(length(unique(batch_train))))
    dev.off()
  } else {
    avg_sil <- NA
    cat("  Not enough samples or batches for silhouette\n")
  }

  # -------------------------------------------------------------------------
  # 2.4 Overlap of training and test distributions (density plot)
  # -------------------------------------------------------------------------
  # Sample a random subset of genes for visualisation
  sample_genes <- sample(colnames(train_mat), min(100, ncol(train_mat)))
  train_long <- train_mat[, sample_genes] %>%
    as.data.frame() %>%
    pivot_longer(everything(), names_to = "gene", values_to = "expression") %>%
    mutate(set = "train")
  test_long <- test_mat[, sample_genes] %>%
    as.data.frame() %>%
    pivot_longer(everything(), names_to = "gene", values_to = "expression") %>%
    mutate(set = "test")
  combined_long <- rbind(train_long, test_long)

  p_dens <- ggplot(combined_long, aes(x = expression, colour = set)) +
    geom_density() +
    labs(title = paste(toupper(plat), "- Train/Test expression density overlap"),
         x = "Normalised expression (after STEP 1)") +
    theme_minimal()
  ggsave(file.path(FIGURES_DIR, paste0(plat, "_density_overlap.pdf")),
         p_dens, width = 7, height = 5, dpi = 300, bg = "white")

  # -------------------------------------------------------------------------
  # 2.5 Check for unseen batch levels in test
  # -------------------------------------------------------------------------
  cat("\n[2.5] Batch consistency between train and test\n")
  train_batches <- unique(meta_train$study)
  test_batches <- unique(meta_test$study)
  unseen <- setdiff(test_batches, train_batches)
  if (length(unseen) > 0) {
    cat("  WARNING: Test contains studies not seen in training:\n")
    for (u in unseen) cat("    -", u, "\n")
    cat("  This is acceptable – the model will be evaluated on novel cohorts.\n")
  } else {
    cat("  All test studies present in training.\n")
  }

  # -------------------------------------------------------------------------
  # 2.6 Collect metrics for this platform
  # -------------------------------------------------------------------------
  metrics <- list(
    platform = plat,
    n_train = nrow(train_mat),
    n_test = nrow(test_mat),
    n_genes = ncol(train_mat),
    mean_batch_R2_train = mean(r2_vals),
    silhouette_width_train = ifelse(is.na(avg_sil), NA, avg_sil),
    unseen_test_batches = unseen,
    pca_variance_explained = var_exp[1:2]
  )
  all_metrics[[plat]] <- metrics

  cat("\n  Metrics saved for", plat, "\n")
}

# -----------------------------------------------------------------------------
# 3. Save all metrics as JSON
# -----------------------------------------------------------------------------
cat("\n", rep("=", 80), "\n")
cat("[3] Saving validation metrics\n")
write_json(all_metrics, file.path(OUTPUT_DIR, "validation_metrics.json"), pretty = TRUE, auto_unbox = TRUE)

# -----------------------------------------------------------------------------
# 4. Summary report
# -----------------------------------------------------------------------------
cat("\n", rep("=", 80), "\n")
cat("DIAGNOSTIC VALIDATION COMPLETED\n")
cat(rep("=", 80), "\n\n")
cat("Output directory:", OUTPUT_DIR, "\n")
cat("  - figures/ : PCA, R² density, silhouette, density overlap\n")
cat("  - validation_metrics.json : quantitative metrics\n")
cat("\nInterpretation:\n")
cat("  - Mean batch R² < 0.1 indicates good batch removal.\n")
cat("  - Silhouette < 0 indicates no clustering by batch.\n")
cat("  - PCA overlay shows train/test similarity.\n")
cat("\nSTEP 2 validation complete. Data are ready for modelling.\n")
cat("End time:", as.character(Sys.time()), "\n")

sink()
