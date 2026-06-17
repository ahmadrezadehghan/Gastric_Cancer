#!/usr/bin/env Rscript
# =============================================================================
# STEP 1: COMPLETE PREPROCESSING – NO LEAKAGE, NATURE MEDICINE READY
# =============================================================================
# - Per‑study stratified train/test split (metadata only)
# - RNA‑seq: voom on training; test: log2(CPM+1e-6) + global quantile harmonisation
# - Microarray: quantile normalisation (training target applied to test)
# - Global quantile harmonisation (target from training only) + optional ComBat+fsva
# - Feature selection by median variance across studies
# - QC plots before/after batch correction, for train and test
# - No imputation before split; NAs handled by gene filtering or column means
# =============================================================================

# -----------------------------------------------------------------------------
# 0. Configuration – EDIT THESE PATHS
# -----------------------------------------------------------------------------
BASE_DIR        <- "E:/GastricCancer-2026"
GENE_SET_PATH   <- file.path(BASE_DIR, "Gene-Set", "Gene_Set.xlsx")
STUDY_INFO_PATH <- file.path(BASE_DIR, "Gene-Set", "study.xlsx")
CLINICAL_PATH   <- file.path(BASE_DIR, "Clinical_Data", "demo.xlsx")
PROBE_MAP_PATH  <- file.path(BASE_DIR, "Annotation", "probe_to_gene.csv")   # optional
OUTPUT_BASE     <- file.path(BASE_DIR, "Prepared_Data")
QC_DIR          <- file.path(OUTPUT_BASE, "QC_plots")
TEMP_DIR        <- file.path(OUTPUT_BASE, "temp")

RANDOM_SEED        <- 42
TRAIN_FRACTION     <- 0.7
MIN_COUNT_RNA      <- 10
MIN_SAMP_FRAC_RNA  <- 0.1
TOP_VARIABLE_GENES <- 20000
BATCH_CORRECT      <- TRUE   # set to FALSE if you rely only on quantile harmonisation

# -----------------------------------------------------------------------------
# 1. Load libraries
# -----------------------------------------------------------------------------
suppressPackageStartupMessages({
  library(readxl)
  library(tidyverse)
  library(DESeq2)
  library(limma)
  library(matrixStats)
  library(jsonlite)
  library(preprocessCore)
  library(sva)
  library(ggplot2)
  library(RColorBrewer)
})

# -----------------------------------------------------------------------------
# 2. Helper functions
# -----------------------------------------------------------------------------

# Collapse probes to genes, resolve duplicate gene symbols (keep max mean)
collapse_to_genes <- function(mat, probe_map = NULL) {
  if (!is.null(probe_map)) {
    common <- intersect(rownames(mat), probe_map$probe_id)
    if (length(common) == 0) stop("No probes in mapping file")
    mat <- mat[common, , drop = FALSE]
    probe_map <- probe_map[match(common, probe_map$probe_id), ]
    genes <- unique(probe_map$gene_symbol)
    new_mat <- matrix(NA, nrow = length(genes), ncol = ncol(mat))
    rownames(new_mat) <- genes
    colnames(new_mat) <- colnames(mat)
    for (gene in genes) {
      probes <- probe_map$probe_id[probe_map$gene_symbol == gene]
      if (length(probes) == 1) {
        new_mat[gene, ] <- mat[probes, ]
      } else {
        means <- rowMeans(mat[probes, , drop = FALSE], na.rm = TRUE)
        best <- probes[which.max(means)]
        new_mat[gene, ] <- mat[best, ]
      }
    }
    mat <- new_mat
  }
  if (any(duplicated(rownames(mat)))) {
    genes <- unique(rownames(mat))
    new_mat <- matrix(NA, nrow = length(genes), ncol = ncol(mat))
    rownames(new_mat) <- genes
    colnames(new_mat) <- colnames(mat)
    for (gene in genes) {
      rows <- which(rownames(mat) == gene)
      if (length(rows) == 1) {
        new_mat[gene, ] <- mat[rows, ]
      } else {
        means <- rowMeans(mat[rows, , drop = FALSE], na.rm = TRUE)
        best <- rows[which.max(means)]
        new_mat[gene, ] <- mat[best, ]
      }
    }
    mat <- new_mat
  }
  return(mat)
}

# QC plots: density + PCA
qc_plots <- function(mat_list, title_prefix, output_dir) {
  all_mat <- do.call(rbind, mat_list)
  study_labels <- rep(names(mat_list), times = sapply(mat_list, nrow))
  dens_data <- data.frame(value = as.vector(all_mat), study = study_labels)
  p_dens <- ggplot(dens_data, aes(x = value, colour = study)) +
    geom_density() + theme_bw() + ggtitle(paste(title_prefix, "- Density"))
  pca <- prcomp(all_mat, center = TRUE, scale. = FALSE)
  pca_df <- data.frame(PC1 = pca$x[,1], PC2 = pca$x[,2], study = study_labels)
  var_exp <- round(100 * summary(pca)$importance[2, 1:2], 1)
  p_pca <- ggplot(pca_df, aes(x = PC1, y = PC2, colour = study)) +
    geom_point(size = 2) + theme_bw() +
    labs(x = paste0("PC1 (", var_exp[1], "%)"), y = paste0("PC2 (", var_exp[2], "%)")) +
    ggtitle(paste(title_prefix, "- PCA"))
  pdf(file.path(output_dir, paste0(gsub(" ", "_", title_prefix), ".pdf")),
      width = 12, height = 5)
  print(p_dens); print(p_pca)
  dev.off()
}

# Check stratified split balance
check_balance <- function(train_ids, test_ids, status_vec) {
  tr_rate <- mean(status_vec[train_ids] == 1, na.rm = TRUE)
  te_rate <- mean(status_vec[test_ids] == 1, na.rm = TRUE)
  cat(sprintf("    Train event rate: %.3f, Test event rate: %.3f, Diff = %.3f\n",
              tr_rate, te_rate, abs(tr_rate - te_rate)))
  if (abs(tr_rate - te_rate) > 0.1) warning("Large difference in event rates")
}

# -----------------------------------------------------------------------------
# 3. Main function for one platform
# -----------------------------------------------------------------------------
run_platform <- function(platform_type) {
  output_dir <- file.path(OUTPUT_BASE, platform_type)
  qc_dir <- file.path(QC_DIR, platform_type)
  temp_dir <- file.path(TEMP_DIR, platform_type)
  for (d in c(output_dir, qc_dir, temp_dir)) dir.create(d, recursive = TRUE, showWarnings = FALSE)

  log_file <- file.path(output_dir, "step1.1_log.txt")
  sink_conn <- file(log_file, open = "wt")
  sink(sink_conn, type = "output", split = TRUE)
  sink(sink_conn, type = "message", append = TRUE)
  on.exit({ sink(type = "output"); sink(type = "message"); close(sink_conn) }, add = TRUE)

  cat(rep("=", 80), "\n")
  cat("STEP 1.1 – PLATFORM:", toupper(platform_type), "\n")
  cat("Start time:", as.character(Sys.time()), "\n")
  cat(rep("=", 80), "\n\n")
  set.seed(RANDOM_SEED)
  saveRDS(.Random.seed, file.path(temp_dir, "random_seed.rds"))

  # -------------------------------------------------------------------------
  # 4. Load expression and study info
  # -------------------------------------------------------------------------
  cat("[1] Loading expression data\n")
  gene_data <- read_excel(GENE_SET_PATH, col_names = TRUE)
  expr_raw <- as.matrix(gene_data[, -1])
  rownames(expr_raw) <- as.character(gene_data[[1]])
  colnames(expr_raw) <- colnames(gene_data)[-1]
  mode(expr_raw) <- "numeric"
  cat("  Raw matrix:", nrow(expr_raw), "genes ×", ncol(expr_raw), "samples\n")

  # Remove all‑zero/all‑NA rows
  keep_rows <- apply(expr_raw, 1, function(x) !all(is.na(x) | x == 0))
  expr_raw <- expr_raw[keep_rows, ]
  cat("  After empty row removal:", nrow(expr_raw), "genes\n")

  # Probe‑to‑gene mapping + duplicate resolution
  if (grepl("^[0-9]+_at", rownames(expr_raw)[1]) || grepl("^[0-9]+$", rownames(expr_raw)[1])) {
    if (!file.exists(PROBE_MAP_PATH)) stop("Probe mapping required but file missing")
    probe_map <- read.csv(PROBE_MAP_PATH, stringsAsFactors = FALSE)
    expr_raw <- collapse_to_genes(expr_raw, probe_map)
    cat("  After probe‑to‑gene mapping:", nrow(expr_raw), "genes\n")
  } else {
    expr_raw <- collapse_to_genes(expr_raw, probe_map = NULL)
    cat("  After duplicate gene resolution:", nrow(expr_raw), "genes\n")
  }

  # Load study mapping
  cat("\n[2] Loading study information\n")
  study_info <- read_excel(STUDY_INFO_PATH)
  study_info$study_name <- study_info$file_name %>%
    str_replace("_series_matrix\\.xlsx$", "") %>%
    str_replace("\\.xlsx$", "") %>%
    str_replace("_series_matrix\\.Asian\\.cohort$", "")

  sample_study <- list()
  for (i in seq_len(nrow(study_info))) {
    study <- study_info$study_name[i]
    samples <- as.character(study_info[i, -1])
    samples <- samples[!is.na(samples) & samples != ""]
    for (s in samples) sample_study[[s]] <- study
  }
  sample_df <- data.frame(sample_id = names(sample_study), study = unlist(sample_study),
                          stringsAsFactors = FALSE)
  sample_df <- sample_df[sample_df$sample_id %in% colnames(expr_raw), ]
  expr_raw <- expr_raw[, sample_df$sample_id, drop = FALSE]
  cat("  Samples with data:", nrow(sample_df), "\n")

  # Filter by platform
  sample_df$platform <- ifelse(grepl("TCGA", sample_df$study, ignore.case = TRUE),
                               "rnaseq", "microarray")
  sample_df <- sample_df[sample_df$platform == platform_type, ]
  if (nrow(sample_df) == 0) { cat("No", platform_type, "studies – skipping\n"); return(NULL) }
  expr_raw <- expr_raw[, sample_df$sample_id, drop = FALSE]
  cat("  Kept", nrow(sample_df), "samples for", platform_type, "\n")

  # -------------------------------------------------------------------------
  # 5. Clinical data for stratification
  # -------------------------------------------------------------------------
  cat("\n[3] Loading clinical data\n")
  survival_status <- NULL; use_strat <- FALSE
  if (file.exists(CLINICAL_PATH)) {
    demo_raw <- read_excel(CLINICAL_PATH, col_names = TRUE)
    demo_mat <- as.matrix(demo_raw[, -1])
    rownames(demo_mat) <- demo_raw[[1]]
    colnames(demo_mat) <- colnames(demo_raw)[-1]
    status_row <- grep("survival status", rownames(demo_mat), ignore.case = TRUE)
    if (length(status_row) == 1) {
      raw_status <- demo_mat[status_row, ]
      survival_status <- setNames(
        ifelse(tolower(raw_status) %in% c("dead", "deceased"), 1,
               ifelse(tolower(raw_status) %in% c("alive", "living"), 0, NA)),
        names(raw_status)
      )
      known_frac <- mean(!is.na(survival_status[sample_df$sample_id]))
      if (known_frac >= 0.7) use_strat <- TRUE
      cat("  Survival known for", round(100*known_frac,1), "%\n")
    }
  }

  # -------------------------------------------------------------------------
  # 6. Train/test split (metadata only)
  # -------------------------------------------------------------------------
  cat("\n[4] Creating train/test split\n")
  train_ids <- c(); test_ids <- c()
  studies <- unique(sample_df$study)
  for (study in studies) {
    study_samples <- sample_df$sample_id[sample_df$study == study]
    n_total <- length(study_samples)
    if (n_total == 1) { train_ids <- c(train_ids, study_samples); next }
    n_train <- max(1, floor(TRAIN_FRACTION * n_total))
    if (n_train == n_total) n_train <- n_total - 1
    if (use_strat) {
      status_vals <- survival_status[study_samples]
      if (all(is.na(status_vals))) {
        train_idx <- sample(seq_len(n_total), n_train, replace = FALSE)
      } else {
        event_idx <- which(status_vals == 1)
        non_idx <- which(status_vals == 0)
        n_event <- floor(TRAIN_FRACTION * length(event_idx))
        n_non <- floor(TRAIN_FRACTION * length(non_idx))
        n_event <- min(n_event, length(event_idx))
        n_non <- min(n_non, length(non_idx))
        train_event <- if (n_event > 0) sample(event_idx, n_event) else integer(0)
        train_non <- if (n_non > 0) sample(non_idx, n_non) else integer(0)
        train_idx <- c(train_event, train_non)
        if (length(train_idx) < n_train) {
          remaining <- setdiff(seq_len(n_total), train_idx)
          train_idx <- c(train_idx, sample(remaining, n_train - length(train_idx), replace = FALSE))
        }
      }
    } else {
      train_idx <- sample(seq_len(n_total), n_train, replace = FALSE)
    }
    train_ids <- c(train_ids, study_samples[train_idx])
    test_ids <- c(test_ids, study_samples[-train_idx])
    if (use_strat) check_balance(train_ids, test_ids, survival_status)
  }
  sample_df$set <- ifelse(sample_df$sample_id %in% train_ids, "train", "test")
  cat("  Total train:", length(train_ids), "samples\n")
  cat("  Total test: ", length(test_ids), "samples\n")
  saveRDS(list(train = train_ids, test = test_ids), file.path(temp_dir, "split_ids.rds"))

  # -------------------------------------------------------------------------
  # 7. Per‑study normalisation (training only, save parameters)
  # -------------------------------------------------------------------------
  cat("\n[5] Normalising training samples\n")
  train_norm_list <- list()
  norm_params <- list()

  for (study in studies) {
    study_train <- sample_df$sample_id[sample_df$study == study & sample_df$set == "train"]
    if (length(study_train) == 0) next
    expr_train <- expr_raw[, study_train, drop = FALSE]  # genes x samples
    cat("  Study:", study, " n_train =", ncol(expr_train), "\n")

    if (platform_type == "rnaseq") {
      # Remove genes with >50% NAs (voom cannot handle NAs)
      na_frac <- apply(expr_train, 1, function(x) mean(is.na(x)))
      expr_train <- expr_train[na_frac <= 0.5, , drop = FALSE]
      # Convert to integer counts (round)
      expr_train <- round(expr_train)
      # Low‑count filter
      min_samp <- max(2, round(MIN_SAMP_FRAC_RNA * ncol(expr_train)))
      keep_genes <- rowSums(expr_train >= MIN_COUNT_RNA, na.rm = TRUE) >= min_samp
      expr_train <- expr_train[keep_genes, , drop = FALSE]
      cat("    Retained", sum(keep_genes), "genes after low‑count filter\n")
      # voom normalisation
      dge <- DGEList(counts = expr_train)
      dge <- calcNormFactors(dge)
      design <- matrix(1, nrow = ncol(expr_train), ncol = 1)  # intercept only
      v <- voom(dge, design = design, plot = FALSE)
      norm_mat <- t(v$E)   # samples x genes
      train_norm_list[[study]] <- norm_mat
      norm_params[[study]] <- list(
        platform = "rnaseq",
        genes_kept = rownames(expr_train),
        voom_trend = v$voom.xy
      )
    } else { # microarray
      # Remove genes with >50% NAs
      na_frac <- apply(expr_train, 1, function(x) mean(is.na(x)))
      expr_train <- expr_train[na_frac <= 0.5, , drop = FALSE]
      # Log2 transform if raw scale
      med_val <- median(expr_train, na.rm = TRUE)
      if (med_val > 10) expr_train <- log2(expr_train + 1)
      # Quantile normalisation (training)
      expr_train_t <- t(expr_train)
      norm_train <- normalizeBetweenArrays(expr_train_t, method = "quantile")
      # Reference quantiles for later test application
      sorted <- apply(norm_train, 2, sort)
      ref_quantiles <- rowMeans(sorted)
      train_norm_list[[study]] <- norm_train
      norm_params[[study]] <- list(
        platform = "microarray",
        already_logged = (med_val <= 10),
        ref_quantiles = ref_quantiles,
        genes_kept = colnames(expr_train)
      )
    }
  }

  # -------------------------------------------------------------------------
  # 8. Apply normalisation to test samples (using saved parameters)
  # -------------------------------------------------------------------------
  cat("\n[6] Applying normalisation to test samples\n")
  test_norm_list <- list()

  for (study in studies) {
    study_test <- sample_df$sample_id[sample_df$study == study & sample_df$set == "test"]
    if (length(study_test) == 0) next
    if (!study %in% names(norm_params)) next
    params <- norm_params[[study]]
    expr_test <- expr_raw[, study_test, drop = FALSE]
    cat("  Study:", study, " n_test =", ncol(expr_test), "\n")

    if (platform_type == "rnaseq") {
      # Subset to genes kept in training
      expr_test <- expr_test[rownames(expr_test) %in% params$genes_kept, , drop = FALSE]
      expr_test <- expr_test[params$genes_kept, , drop = FALSE]
      # Compute log2(CPM + 1e-6) – safe, no -Inf, no leakage
      lib_sizes <- colSums(expr_test, na.rm = TRUE)
      cpm <- sweep(expr_test, 2, lib_sizes / 1e6, FUN = "/")
      logcpm <- log2(cpm + 1e-6)   # 1e-6 avoids -Inf
      # Transpose to samples x genes
      test_norm_list[[study]] <- t(logcpm)
    } else { # microarray
      if (!params$already_logged) {
        med_val <- median(expr_test, na.rm = TRUE)
        if (med_val > 10) expr_test <- log2(expr_test + 1)
      }
      # Subset to training genes
      expr_test <- expr_test[rownames(expr_test) %in% params$genes_kept, , drop = FALSE]
      expr_test <- expr_test[params$genes_kept, , drop = FALSE]
      expr_test_t <- t(expr_test)
      norm_test <- normalize.quantiles.use.target(expr_test_t, target = params$ref_quantiles)
      test_norm_list[[study]] <- norm_test
    }
  }

  if (length(test_norm_list) == 0) {
    cat("WARNING: No test samples after normalisation. Test output will be empty.\n")
  }

  # -------------------------------------------------------------------------
  # 9. Global quantile harmonisation (target from training only)
  # -------------------------------------------------------------------------
  cat("\n[7] Global quantile harmonisation (training target applied to test)\n")
  train_combined <- do.call(rbind, train_norm_list)
  # Compute reference quantiles from training
  train_sorted <- apply(train_combined, 2, sort)
  global_target <- rowMeans(train_sorted)
  # Apply to training (in place)
  train_harm <- normalize.quantiles.use.target(train_combined, target = global_target)
  rownames(train_harm) <- rownames(train_combined)
  colnames(train_harm) <- colnames(train_combined)
  # Re‑split training list for QC
  train_harm_list <- list()
  start <- 1
  for (nm in names(train_norm_list)) {
    n <- nrow(train_norm_list[[nm]])
    train_harm_list[[nm]] <- train_harm[start:(start+n-1), , drop = FALSE]
    start <- start + n
  }
  # Apply to test
  if (length(test_norm_list) > 0) {
    test_combined <- do.call(rbind, test_norm_list)
    test_harm <- normalize.quantiles.use.target(test_combined, target = global_target)
    rownames(test_harm) <- rownames(test_combined)
    colnames(test_harm) <- colnames(test_combined)
    test_harm_list <- list()
    start <- 1
    for (nm in names(test_norm_list)) {
      n <- nrow(test_norm_list[[nm]])
      test_harm_list[[nm]] <- test_harm[start:(start+n-1), , drop = FALSE]
      start <- start + n
    }
  } else {
    test_harm_list <- list()
  }

  # QC before batch correction
  qc_plots(train_harm_list, paste0(platform_type, "_after_harmonisation"), qc_dir)
  if (length(test_harm_list) > 0)
    qc_plots(test_harm_list, paste0(platform_type, "_test_after_harmonisation"), qc_dir)

  # -------------------------------------------------------------------------
  # 10. Batch correction (ComBat on training, fsva to test)
  # -------------------------------------------------------------------------
  cat("\n[8] Batch correction across studies\n")
  if (BATCH_CORRECT && length(unique(sample_df$study[sample_df$set == "train"])) > 1) {
    train_mat <- do.call(rbind, train_harm_list)
    batch_train <- sample_df$study[match(rownames(train_mat), sample_df$sample_id)]
    # ComBat (genes x samples)
    train_combat <- ComBat(t(train_mat), batch = batch_train, mod = NULL, par.prior = TRUE)
    train_corrected <- t(train_combat)
    # Apply to test via fsva
    if (length(test_harm_list) > 0) {
      test_mat <- do.call(rbind, test_harm_list)
      fsva_obj <- fsva(dbdat = t(train_mat), mod = NULL, sv = NULL, newdat = t(test_mat))
      test_corrected <- t(fsva_obj$new)
    } else {
      test_corrected <- matrix(0, nrow = 0, ncol = ncol(train_corrected))
    }
    cat("  Batch correction applied (ComBat + fsva)\n")
  } else {
    train_corrected <- do.call(rbind, train_harm_list)
    test_corrected <- if (length(test_harm_list) > 0) do.call(rbind, test_harm_list) else matrix(0, nrow = 0, ncol = ncol(train_corrected))
  }

  # QC after batch correction
  qc_plots(list(train_corrected = train_corrected), paste0(platform_type, "_train_after_batch"), qc_dir)
  if (nrow(test_corrected) > 0)
    qc_plots(list(test_corrected = test_corrected), paste0(platform_type, "_test_after_batch"), qc_dir)

  # -------------------------------------------------------------------------
  # 11. Feature selection: median variance across studies
  # -------------------------------------------------------------------------
  cat("\n[9] Feature selection (median variance across studies)\n")
  studies_train <- unique(sample_df$study[sample_df$set == "train"])
  gene_var_matrix <- matrix(NA, nrow = ncol(train_corrected), ncol = length(studies_train))
  rownames(gene_var_matrix) <- colnames(train_corrected)
  colnames(gene_var_matrix) <- studies_train
  for (st in studies_train) {
    samples_st <- sample_df$sample_id[sample_df$study == st & sample_df$set == "train"]
    if (length(samples_st) > 1) {
      mat_st <- train_corrected[rownames(train_corrected) %in% samples_st, , drop = FALSE]
      gene_var_matrix[, st] <- colVars(mat_st, na.rm = TRUE)
    }
  }
  median_var <- apply(gene_var_matrix, 1, median, na.rm = TRUE)
  median_var <- median_var[!is.na(median_var) & median_var > 0]
  if (!is.null(TOP_VARIABLE_GENES) && length(median_var) > TOP_VARIABLE_GENES) {
    top_genes <- names(sort(median_var, decreasing = TRUE)[1:TOP_VARIABLE_GENES])
  } else {
    top_genes <- names(median_var)
  }
  train_final <- train_corrected[, top_genes, drop = FALSE]
  cat("  Final training features:", ncol(train_final), "\n")

  # -------------------------------------------------------------------------
  # 12. Align test data and impute missing with training column means
  # -------------------------------------------------------------------------
  cat("\n[10] Aligning test data\n")
  if (nrow(test_corrected) > 0) {
    common <- intersect(colnames(test_corrected), colnames(train_final))
    if (length(common) == 0) stop("No common genes between train and test")
    test_sub <- test_corrected[, common, drop = FALSE]
    missing <- setdiff(colnames(train_final), colnames(test_sub))
    if (length(missing) > 0) {
      na_mat <- matrix(NA, nrow = nrow(test_sub), ncol = length(missing))
      colnames(na_mat) <- missing
      test_sub <- cbind(test_sub, na_mat)
    }
    test_sub <- test_sub[, colnames(train_final), drop = FALSE]
    train_means <- colMeans(train_final, na.rm = TRUE)
    for (j in 1:ncol(test_sub)) {
      na_idx <- is.na(test_sub[, j])
      if (any(na_idx)) test_sub[na_idx, j] <- train_means[j]
    }
    test_final <- test_sub
  } else {
    test_final <- matrix(0, nrow = 0, ncol = ncol(train_final))
  }
  cat("  Final test set:", nrow(test_final), "samples ×", ncol(test_final), "genes\n")

  # -------------------------------------------------------------------------
  # 13. Save outputs
  # -------------------------------------------------------------------------
  cat("\n[11] Saving outputs\n")
  write.csv(as.data.frame(train_final), file.path(output_dir, "train_data.csv"), row.names = TRUE)
  write.csv(as.data.frame(test_final),  file.path(output_dir, "test_data.csv"),  row.names = TRUE)

  sample_metadata <- sample_df[sample_df$sample_id %in% c(rownames(train_final), rownames(test_final)), ]
  write.csv(sample_metadata, file.path(output_dir, "sample_metadata.csv"), row.names = FALSE)

  study_platform <- unique(sample_metadata[, c("study", "platform")])
  write.csv(study_platform, file.path(output_dir, "study_platform.csv"), row.names = FALSE)

  # Save parameters for reproducibility
  saveRDS(norm_params, file.path(temp_dir, "norm_params.rds"))
  if (BATCH_CORRECT && exists("fsva_obj")) saveRDS(fsva_obj, file.path(temp_dir, "fsva_obj.rds"))
  saveRDS(train_means, file.path(temp_dir, "train_means.rds"))
  saveRDS(top_genes, file.path(temp_dir, "top_genes.rds"))
  saveRDS(global_target, file.path(temp_dir, "quantile_target.rds"))

  metadata <- list(
    pipeline = "Step 1.1 – no leakage, complete",
    platform = platform_type,
    date = as.character(Sys.time()),
    split_method = ifelse(use_strat, "stratified by survival", "random"),
    train_fraction = TRAIN_FRACTION,
    random_seed = RANDOM_SEED,
    n_train = nrow(train_final),
    n_test = nrow(test_final),
    n_genes = ncol(train_final),
    batch_corrected = BATCH_CORRECT,
    top_variable_genes = ifelse(is.null(TOP_VARIABLE_GENES), "all", TOP_VARIABLE_GENES),
    warning = "Internal train/test split does NOT replace external validation. For Nature Medicine, an independent external cohort is required."
  )
  write_json(metadata, file.path(output_dir, "metadata_step1.1.json"), pretty = TRUE, auto_unbox = TRUE)
  writeLines(capture.output(sessionInfo()), file.path(output_dir, "sessionInfo.txt"))

  cat("\n", rep("=", 80), "\n")
  cat("STEP 1.1 COMPLETED FOR PLATFORM:", toupper(platform_type), "\n")
  cat("End time:", as.character(Sys.time()), "\n")
  return(metadata)
}

# -----------------------------------------------------------------------------
# 14. Execute for each platform present
# -----------------------------------------------------------------------------
# Quick platform detection using first few rows
gene_data_tmp <- read_excel(GENE_SET_PATH, n_max = 1)
samples_tmp <- colnames(gene_data_tmp)[-1]
study_info_tmp <- read_excel(STUDY_INFO_PATH)
study_info_tmp$study_name <- study_info_tmp$file_name %>%
  str_replace("_series_matrix\\.xlsx$", "") %>% str_replace("\\.xlsx$", "")
sample_study_tmp <- unlist(lapply(seq_len(nrow(study_info_tmp)), function(i) {
  rep(study_info_tmp$study_name[i], sum(!is.na(study_info_tmp[i, -1]) & study_info_tmp[i, -1] != ""))
}))
platforms_present <- unique(ifelse(grepl("TCGA", sample_study_tmp, ignore.case = TRUE), "rnaseq", "microarray"))

for (p in platforms_present) {
  tryCatch(run_platform(p), error = function(e) {
    cat("ERROR in platform", p, ":", e$message, "\n")
    sink(type = "output"); sink(type = "message")
  })
}

cat("\nAll done.\n")
cat("IMPORTANT: Internal train/test split is NOT external validation.\n")
cat("To meet Nature Medicine standards, validate on an independent cohort.\n")
