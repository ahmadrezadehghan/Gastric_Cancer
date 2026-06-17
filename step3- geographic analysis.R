#!/usr/bin/env Rscript
# =============================================================================
# STEP 3: COUNTRY & REGION‑SPECIFIC SURVIVAL ANALYSIS (FIXED, NO DATA LEAKAGE)
# =============================================================================
# This script performs two integrated analyses:
#
#   PART A: Country effects on survival
#           - World maps of sample distribution (binary + gradient)
#           - Country‑level hazard ratios (Cox, reference = largest country)
#           - Subregion Kaplan‑Meier curves and Cox regression
#
#   PART B: Region‑specific prognostic gene discovery (TRAINING SET ONLY)
#           - Univariate Cox per subregion (top 2000 variable genes)
#           - Volcano plots with top gene labels
#           - Overlap analysis (Upset plot, bar chart)
#           - Pairwise HR scatter plots
#           - Hallmark pathway enrichment (hypergeometric test, using msigdbr)
#           - Heatmap of common prognostic genes
#           - Summary table of all p‑values
#
# All outputs (CSV, PNG, PDF) are saved to dedicated directories.
# =============================================================================

# -----------------------------------------------------------------------------
# 0. User configuration (EDIT THESE PATHS)
# -----------------------------------------------------------------------------
BASE_DIR <- "E:/GastricCancer-2026"
PLATFORM <- "rnaseq"   # or "microarray" – must match Step 1.1 output subfolder
CLINICAL_FILE <- file.path(BASE_DIR, "Clinical_Data", "Demographical_Gastric_UsedRows_Only.xlsx")
EXPR_TRAIN_FILE <- file.path(BASE_DIR, "Prepared_Data", PLATFORM, "train_data.csv")
MAPPING_FILE <- file.path(BASE_DIR, "Gene-Set", "frame.xlsx")
OUTPUT_DIR <- file.path(BASE_DIR, "Country_Region_Analysis")
FIG_DIR <- file.path(OUTPUT_DIR, "figures")

# Analysis parameters
MIN_SAMPLES_PER_REGION <- 30      # minimum samples for region‑specific Cox
MIN_EVENTS_PER_REGION <- 10       # minimum events for region‑specific Cox
TOP_VARIABLE_GENES <- 2000        # number of most variable genes to test per region
N_TOP_GENES_VOLCANO <- 10         # number of top genes to label on volcano plots
FDR_THRESHOLD <- 0.05             # significance threshold for FDR

# -----------------------------------------------------------------------------
# 1. Load required packages (install if missing)
# -----------------------------------------------------------------------------
required_pkgs <- c(
  "survival", "survminer", "dplyr", "ggplot2", "readxl", "tidyr", "stringr",
  "forcats", "RColorBrewer", "rnaturalearth", "rnaturalearthdata", "sf",
  "viridis", "pheatmap", "UpSetR", "ggrepel", "matrixStats", "jsonlite"
)
# msigdbr is optional but recommended for Hallmark sets
if (!require("msigdbr", quietly = TRUE)) {
  cat("msigdbr not installed. Will use built‑in fallback Hallmark sets.\n")
} else {
  required_pkgs <- c(required_pkgs, "msigdbr")
}

for (pkg in required_pkgs) {
  if (!require(pkg, character.only = TRUE)) {
    install.packages(pkg, dependencies = TRUE)
    library(pkg, character.only = TRUE)
  }
}

# Create output directories
if (!dir.exists(OUTPUT_DIR)) dir.create(OUTPUT_DIR, recursive = TRUE)
if (!dir.exists(FIG_DIR)) dir.create(FIG_DIR, recursive = TRUE)

# Start logging
log_file <- file.path(OUTPUT_DIR, "analysis_log.txt")
sink(log_file, append = FALSE, split = TRUE)

cat(rep("=", 80), "\n")
cat("STEP 3: COUNTRY & REGION‑SPECIFIC SURVIVAL ANALYSIS\n")
cat("Start time:", date(), "\n")
cat(rep("=", 80), "\n\n")

# Helper function to save plots as both PDF and PNG
save_plot <- function(plot, filename, width = 12, height = 6, dpi = 300) {
  ggsave(paste0(filename, ".pdf"), plot = plot, width = width, height = height, dpi = dpi)
  ggsave(paste0(filename, ".png"), plot = plot, width = width, height = height, dpi = dpi, bg = "white")
  cat("  Saved:", filename, "(PDF and PNG)\n")
}

# Helper to normalise sample IDs (TCGA uses dashes, others dots)
normalize_id <- function(x) {
  x <- as.character(x)
  ifelse(grepl("TCGA", x), gsub("\\.", "-", x), x)
}

# -----------------------------------------------------------------------------
# PART A: COUNTRY EFFECTS ON SURVIVAL (using clinical data only)
# -----------------------------------------------------------------------------
cat("\n", rep("=", 80), "\n")
cat("PART A: COUNTRY EFFECTS ON SURVIVAL\n")
cat(rep("=", 80), "\n\n")

# A1. Load clinical data
cat("[A1] Loading clinical data from", CLINICAL_FILE, "\n")
if (!file.exists(CLINICAL_FILE)) stop("Clinical file not found: ", CLINICAL_FILE)

clin_raw <- read_excel(CLINICAL_FILE, col_names = TRUE)
feature_names <- clin_raw[[1]]
clin_matrix <- as.matrix(clin_raw[, -1])
rownames(clin_matrix) <- feature_names
colnames(clin_matrix) <- colnames(clin_raw)[-1]

clinical_df <- data.frame(
  sample_id = colnames(clin_matrix),
  OS_time = as.numeric(clin_matrix["Overall_survival(Months)", ]),
  OS_status = ifelse(tolower(clin_matrix["Survival status", ]) %in% c("dead", "deceased"), 1, 0),
  Country = clin_matrix["Country", ],
  stringsAsFactors = FALSE
)

clinical_df <- clinical_df[!is.na(clinical_df$OS_time) & !is.na(clinical_df$OS_status) & clinical_df$OS_time > 0, ]
cat("  Samples with valid survival:", nrow(clinical_df), "\n")

# A2. Clean country names
cat("\n[A2] Cleaning country names...\n")
clean_country <- function(x) {
  if (is.na(x)) return(NA)
  if (grepl(",", x)) {
    parts <- strsplit(x, ",")[[1]]
    x <- trimws(parts[length(parts)])
  } else {
    x <- trimws(x)
  }
  x <- tolower(x)
  x <- gsub("south korea|seoul|goyang|daegu", "south korea", x)
  x <- gsub("usa|united states|united states of america", "usa", x)
  x <- gsub("uk|united kingdom", "uk", x)
  x <- gsub("germany|deutschland", "germany", x)
  x <- gsub("china|beijing|shanghai", "china", x)
  x <- gsub("japan|tokyo", "japan", x)
  x <- gsub("australia|melbourne", "australia", x)
  x <- gsub("brazil|são paulo", "brazil", x)
  x <- gsub("mexico|ciudad de méxico", "mexico", x)
  x <- gsub("switzerland|basel", "switzerland", x)
  x <- gsub("france|paris", "france", x)
  x <- gsub("netherlands", "netherlands", x)
  x <- gsub("austria", "austria", x)
  x <- gsub("belgium", "belgium", x)
  x <- gsub("hungary", "hungary", x)
  x <- gsub("romania", "romania", x)
  x <- gsub("singapore", "singapore", x)
  tools::toTitleCase(x)
}
clinical_df$Country_clean <- sapply(clinical_df$Country, clean_country)
cat("  Unique countries after cleaning:\n")
print(unique(clinical_df$Country_clean))

# A3. World map – sample counts (with error handling for rnaturalearth)
cat("\n[A3] Generating world maps for sample distribution...\n")
country_data <- as.data.frame(table(clinical_df$Country_clean), stringsAsFactors = FALSE)
colnames(country_data) <- c("country", "n_samples")
country_data$country <- as.character(country_data$country)

world <- NULL
map_available <- FALSE
tryCatch({
  world <- ne_countries(scale = "medium", returnclass = "sf")
  map_available <- TRUE
}, error = function(e) {
  cat("  Warning: Could not download world map from rnaturalearth. Skipping maps.\n")
})

if (map_available) {
  name_fix <- c("Usa" = "United States of America", "South Korea" = "South Korea",
                "China" = "China", "Brazil" = "Brazil", "Switzerland" = "Switzerland")
  country_data$world_name <- country_data$country
  for (i in 1:nrow(country_data)) {
    cnt <- country_data$country[i]
    if (cnt %in% names(name_fix)) country_data$world_name[i] <- name_fix[cnt]
  }
  world_merged <- merge(world, country_data, by.x = "admin", by.y = "world_name", all.x = TRUE)
  world_merged$n_samples[is.na(world_merged$n_samples)] <- NA
  world_merged$has_samples <- ifelse(!is.na(world_merged$n_samples) & world_merged$n_samples > 0, TRUE, NA)

  map_binary <- ggplot(world_merged) +
    geom_sf(aes(fill = has_samples), color = "gray90", size = 0.1) +
    scale_fill_manual(values = c("TRUE" = "steelblue"), na.value = "white", name = "Has samples") +
    labs(title = "Geographic distribution of samples", subtitle = "Blue = at least one sample") +
    theme_minimal() + theme(plot.title = element_text(hjust = 0.5, face = "bold"))
  save_plot(map_binary, file.path(OUTPUT_DIR, "World_map_sample_counts"))

  map_gradient <- ggplot(world_merged) +
    geom_sf(aes(fill = n_samples), color = "gray90", size = 0.1) +
    scale_fill_gradient(low = "lightblue", high = "darkblue", na.value = "white",
                        trans = "log1p", name = "Number of samples\n(log scale)") +
    labs(title = "Geographic distribution (gradient)") +
    theme_minimal() + theme(plot.title = element_text(hjust = 0.5, face = "bold"))
  save_plot(map_gradient, file.path(OUTPUT_DIR, "World_map_sample_counts_gradient"))
} else {
  cat("  Maps skipped due to missing rnaturalearth data.\n")
}

# A4. Country‑level hazard ratios (countries with ≥5 samples)
cat("\n[A4] Computing country‑level hazard ratios...\n")
country_counts <- country_data[country_data$n_samples >= 5, ]
if (nrow(country_counts) > 1 && map_available) {
  clinical_df$Country_factor <- factor(clinical_df$Country_clean)
  ref_country <- as.character(country_counts$country[which.max(country_counts$n_samples)])
  clinical_df$Country_factor <- relevel(clinical_df$Country_factor, ref = ref_country)
  cox_country <- coxph(Surv(OS_time, OS_status) ~ Country_factor, data = clinical_df)
  coef_sum <- summary(cox_country)$coefficients
  if (nrow(coef_sum) > 0) {
    hr_df <- data.frame(
      country = gsub("Country_factor", "", rownames(coef_sum)),
      HR = exp(coef_sum[, "coef"]),
      lower_95 = exp(coef_sum[, "coef"] - 1.96 * coef_sum[, "se(coef)"]),
      upper_95 = exp(coef_sum[, "coef"] + 1.96 * coef_sum[, "se(coef)"]),
      p = coef_sum[, "Pr(>|z|)"],
      stringsAsFactors = FALSE
    )
    hr_df$world_name <- hr_df$country
    for (i in 1:nrow(hr_df)) {
      if (hr_df$country[i] %in% names(name_fix))
        hr_df$world_name[i] <- name_fix[hr_df$country[i]]
    }
    world_hr <- merge(world, hr_df, by.x = "admin", by.y = "world_name", all.x = TRUE)
    world_hr$HR[is.na(world_hr$HR)] <- NA
    world_hr$HR[world_hr$admin == name_fix[ref_country]] <- NA
    map_hr <- ggplot(world_hr) +
      geom_sf(aes(fill = HR), color = "gray90", size = 0.1) +
      scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 1,
                           na.value = "white", name = "Hazard Ratio") +
      labs(title = "Hazard Ratio by country", subtitle = paste("Reference:", ref_country)) +
      theme_minimal() + theme(plot.title = element_text(hjust = 0.5, face = "bold"))
    save_plot(map_hr, file.path(OUTPUT_DIR, "World_map_HR"))
    write.csv(hr_df, file.path(OUTPUT_DIR, "Country_HR_results.csv"), row.names = FALSE)
  }
} else {
  cat("  Not enough countries with ≥5 samples for HR map.\n")
}

# A5. Subregion analysis (manual mapping – expanded)
cat("\n[A5] Subregion analysis (manual grouping)...\n")
country_to_subregion <- list(
  "Usa" = "Northern America", "Canada" = "Northern America",
  "Mexico" = "Central America", "Brazil" = "South America",
  "Argentina" = "South America", "Chile" = "South America",
  "China" = "Eastern Asia", "South Korea" = "Eastern Asia", "Japan" = "Eastern Asia",
  "Switzerland" = "Western Europe", "France" = "Western Europe", "Germany" = "Western Europe",
  "Austria" = "Western Europe", "Belgium" = "Western Europe", "Netherlands" = "Western Europe",
  "United Kingdom" = "Northern Europe", "Australia" = "Oceania", "Singapore" = "Southeast Asia"
)
clinical_df$Subregion <- sapply(clinical_df$Country_clean, function(cnt) {
  if (is.na(cnt)) return(NA)
  reg <- country_to_subregion[[cnt]]
  if (is.null(reg)) NA else reg
})
clinical_df <- clinical_df[!is.na(clinical_df$Subregion), ]
cat("  Subregion distribution:\n"); print(table(clinical_df$Subregion))

# Kaplan‑Meier by subregion
if (length(unique(clinical_df$Subregion)) > 1) {
  clinical_df$Subregion <- factor(clinical_df$Subregion)
  ref_sub <- names(sort(table(clinical_df$Subregion), decreasing = TRUE))[1]
  clinical_df$Subregion <- relevel(clinical_df$Subregion, ref = ref_sub)
  fit_km <- survfit(Surv(OS_time, OS_status) ~ Subregion, data = clinical_df)
  km_plot <- ggsurvplot(fit_km, data = clinical_df, pval = TRUE, conf.int = TRUE,
                        risk.table = TRUE, xlab = "Months",
                        title = paste("Overall Survival by Subregion (ref =", ref_sub, ")"),
                        legend.title = "Subregion", palette = "Set1")
  save_plot(km_plot$plot, file.path(OUTPUT_DIR, "KM_by_Subregion"), width = 10, height = 7)

  # Cox regression for subregion
  cox_sub <- coxph(Surv(OS_time, OS_status) ~ Subregion, data = clinical_df)
  cox_sub_res <- data.frame(
    Subregion = names(coef(cox_sub)),
    HR = exp(coef(cox_sub)),
    lower_95 = exp(confint(cox_sub)[,1]),
    upper_95 = exp(confint(cox_sub)[,2]),
    p = coef(summary(cox_sub))[, "Pr(>|z|)"]
  )
  write.csv(cox_sub_res, file.path(OUTPUT_DIR, "Cox_subregion_results.csv"), row.names = FALSE)

  # Optional: check proportional hazards assumption (commented)
  # ph_test <- cox.zph(cox_sub)
  # print(ph_test)
  # if (ph_test$table[,"p"] < 0.05) cat("  Warning: Proportional hazards assumption violated for some subregions.\n")
} else {
  cat("  Only one subregion; skipping subregion survival analysis.\n")
}

# Bar plot of subregion sample counts
bar_data <- as.data.frame(table(clinical_df$Subregion))
colnames(bar_data) <- c("Subregion", "Count")
bar_plot <- ggplot(bar_data, aes(x = reorder(Subregion, -Count), y = Count, fill = Subregion)) +
  geom_bar(stat = "identity") + scale_fill_brewer(palette = "Set2") +
  labs(title = "Sample distribution by subregion", x = "", y = "Number of samples") +
  theme_minimal() + theme(axis.text.x = element_text(angle = 45, hjust = 1))
save_plot(bar_plot, file.path(OUTPUT_DIR, "Subregion_distribution"), width = 8, height = 5)

# Save cleaned clinical data for Part B
clin_clean <- clinical_df[, c("sample_id", "OS_time", "OS_status", "Subregion")]
rownames(clin_clean) <- clin_clean$sample_id
write.csv(clin_clean, file.path(OUTPUT_DIR, "clin_clean_with_subregion.csv"), row.names = FALSE)

# -----------------------------------------------------------------------------
# PART B: REGION‑SPECIFIC PROGNOSTIC GENE DISCOVERY (TRAINING SET ONLY)
# -----------------------------------------------------------------------------
cat("\n", rep("=", 80), "\n")
cat("PART B: REGION‑SPECIFIC PROGNOSTIC GENE DISCOVERY (TRAINING SET ONLY)\n")
cat(rep("=", 80), "\n\n")

# B1. Load training expression data (batch‑corrected)
cat("[B1] Loading training expression data from", EXPR_TRAIN_FILE, "\n")
if (!file.exists(EXPR_TRAIN_FILE)) stop("Train expression file not found. Run Step 1.1 first.")
train_expr <- read.csv(EXPR_TRAIN_FILE, row.names = 1, check.names = FALSE)
rownames(train_expr) <- normalize_id(rownames(train_expr))
cat("  Training expression matrix:", nrow(train_expr), "samples ×", ncol(train_expr), "genes\n")

# B2. Load clinical subregion data (from Part A) and match samples
cat("\n[B2] Matching samples with subregion information...\n")
clin_sub <- read.csv(file.path(OUTPUT_DIR, "clin_clean_with_subregion.csv"), row.names = 1)
rownames(clin_sub) <- normalize_id(rownames(clin_sub))
common <- intersect(rownames(train_expr), rownames(clin_sub))
cat("  Common samples:", length(common), "\n")
if (length(common) == 0) stop("No overlapping samples between expression and clinical data.")
train_expr <- train_expr[common, , drop = FALSE]
clin_sub <- clin_sub[common, ]

# Keep only subregions with sufficient samples/events
keep_regions <- names(which(table(clin_sub$Subregion) >= MIN_SAMPLES_PER_REGION))
keep_regions <- keep_regions[sapply(keep_regions, function(r) {
  sum(clin_sub$OS_status[clin_sub$Subregion == r]) >= MIN_EVENTS_PER_REGION
})]
cat("  Regions meeting criteria (≥", MIN_SAMPLES_PER_REGION, " samples, ≥", MIN_EVENTS_PER_REGION, " events):\n")
print(keep_regions)
if (length(keep_regions) == 0) stop("No region meets the minimum sample/event criteria.")
clin_sub <- clin_sub[clin_sub$Subregion %in% keep_regions, , drop = FALSE]
train_expr <- train_expr[rownames(clin_sub), , drop = FALSE]
clin_sub$Subregion <- droplevels(clin_sub$Subregion)

# B3. Load gene symbol mapping
cat("\n[B3] Loading gene symbol mapping from", MAPPING_FILE, "\n")
if (!file.exists(MAPPING_FILE)) stop("Mapping file not found: ", MAPPING_FILE)
mapping_df <- read_excel(MAPPING_FILE, col_names = FALSE)
colnames(mapping_df) <- c("Symbol", "ID")
mapping_df <- mapping_df[!is.na(mapping_df$Symbol) & !is.na(mapping_df$ID), ]
id_to_symbol <- setNames(mapping_df$Symbol, mapping_df$ID)
symbol_to_id <- setNames(mapping_df$ID, mapping_df$Symbol)
map_gene_id <- function(ids) {
  sapply(ids, function(id) if (id %in% names(id_to_symbol)) id_to_symbol[[id]] else id)
}
cat("  Loaded", length(id_to_symbol), "mappings\n")

# B4. Univariate Cox per region (top variable genes)
cat("\n[B4] Running univariate Cox per region (top", TOP_VARIABLE_GENES, "variable genes)...\n")
region_results <- list()

for (reg in levels(clin_sub$Subregion)) {
  cat("\n=== Region:", reg, "===\n")
  keep <- clin_sub$Subregion == reg
  expr_reg <- as.matrix(train_expr[keep, , drop = FALSE])
  clin_reg <- clin_sub[keep, ]
  cat("  n =", nrow(clin_reg), ", events =", sum(clin_reg$OS_status), "\n")

  # Remove zero‑variance genes
  gene_var <- matrixStats::colVars(expr_reg, na.rm = TRUE)
  valid <- which(!is.na(gene_var) & gene_var > 0)
  expr_reg <- expr_reg[, valid, drop = FALSE]
  gene_var <- gene_var[valid]

  # Select top variable genes
  n_genes <- min(TOP_VARIABLE_GENES, length(gene_var))
  top_idx <- order(gene_var, decreasing = TRUE)[1:n_genes]
  top_ids <- names(gene_var)[top_idx]
  expr_sub <- expr_reg[, top_ids, drop = FALSE]

  # Univariate Cox
  res <- data.frame(Gene_id = top_ids, HR = NA, lower = NA, upper = NA, p = NA, stringsAsFactors = FALSE)
  for (i in 1:ncol(expr_sub)) {
    cox <- tryCatch(coxph(Surv(clin_reg$OS_time, clin_reg$OS_status) ~ expr_sub[, i]),
                    error = function(e) NULL)
    if (!is.null(cox)) {
      res$HR[i] <- exp(coef(cox))
      ci <- exp(confint(cox))
      res$lower[i] <- ci[1]
      res$upper[i] <- ci[2]
      res$p[i] <- summary(cox)$coefficients[, "Pr(>|z|)"]
    }
  }
  res$FDR <- p.adjust(res$p, method = "BH")
  res <- res[order(res$p), ]
  res$Gene <- map_gene_id(res$Gene_id)
  res_out <- res[, c("Gene", "Gene_id", "HR", "lower", "upper", "p", "FDR")]

  write.csv(res_out, file.path(OUTPUT_DIR, paste0(reg, "_cox_results.csv")), row.names = FALSE)
  sig_genes <- res_out$Gene[res_out$FDR < FDR_THRESHOLD & !is.na(res_out$FDR)]
  cat("  Significant genes (FDR <", FDR_THRESHOLD, "):", length(sig_genes), "\n")

  region_results[[reg]] <- list(res = res_out, sig_genes = sig_genes, expr = expr_sub, clin = clin_reg)
}

if (length(region_results) == 0) stop("No results generated for any region.")

# B5. Volcano plots for each region
cat("\n[B5] Generating volcano plots...\n")
for (reg in names(region_results)) {
  res <- region_results[[reg]]$res
  res$logP <- -log10(res$p)
  res$logHR <- log2(res$HR)
  res$Significant <- ifelse(res$FDR < FDR_THRESHOLD, "FDR<0.05", "Not significant")
  top_genes <- head(res[order(res$p), "Gene"], N_TOP_GENES_VOLCANO)
  res$Label <- ifelse(res$Gene %in% top_genes, res$Gene, NA)

  p <- ggplot(res, aes(x = logHR, y = logP, color = Significant, label = Label)) +
    geom_point(alpha = 0.6, size = 1.5) +
    scale_color_manual(values = c("FDR<0.05" = "red", "Not significant" = "gray70")) +
    geom_vline(xintercept = 0, linetype = "dashed") +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
    geom_text_repel(size = 3, max.overlaps = 15, na.rm = TRUE) +
    labs(x = expression(log[2] ~ "Hazard Ratio"), y = expression(-log[10] ~ "P-value"),
         title = paste(reg, "- Gene expression vs survival"),
         subtitle = paste(sum(res$FDR < FDR_THRESHOLD, na.rm = TRUE), "significant genes")) +
    theme_classic(base_size = 12) + theme(legend.position = "bottom")
  ggsave(file.path(FIG_DIR, paste0(gsub(" ", "_", reg), "_volcano.png")), p, width = 8, height = 6, dpi = 300)
}

# B6. Overlap analysis (significant genes across regions)
cat("\n[B6] Overlap analysis of significant genes...\n")
gene_lists <- lapply(region_results, function(x) x$sig_genes)
gene_lists <- gene_lists[sapply(gene_lists, length) > 0]
if (length(gene_lists) >= 2) {
  reg_names <- names(gene_lists)
  overlap_mat <- matrix(0, nrow = length(reg_names), ncol = length(reg_names), dimnames = list(reg_names, reg_names))
  for (i in 1:length(reg_names)) {
    for (j in 1:length(reg_names)) {
      overlap_mat[i, j] <- length(intersect(gene_lists[[reg_names[i]]], gene_lists[[reg_names[j]]]))
    }
  }
  write.csv(overlap_mat, file.path(OUTPUT_DIR, "gene_overlap_matrix.csv"), row.names = TRUE)

  # Upset plot
  if (require(UpSetR, quietly = TRUE)) {
    upset_data <- list()
    for (reg in reg_names) upset_data[[reg]] <- gene_lists[[reg]]
    png(file.path(FIG_DIR, "upset_plot.png"), width = 10, height = 6, units = "in", res = 300)
    upset(fromList(upset_data), order.by = "freq", nsets = length(reg_names),
          mainbar.y.label = "Gene set intersection size")
    dev.off()
  }

  # Bar plot: unique vs shared
  all_genes <- unique(unlist(gene_lists))
  common_all <- Reduce(intersect, gene_lists)
  only_one <- sapply(gene_lists, function(x) length(setdiff(x, common_all)))
  bar_df <- data.frame(Region = names(only_one), Unique = only_one, Shared = length(common_all))
  bar_long <- pivot_longer(bar_df, cols = c(Unique, Shared), names_to = "Type", values_to = "Count")
  p_bar <- ggplot(bar_long, aes(x = Region, y = Count, fill = Type)) +
    geom_bar(stat = "identity", position = "dodge") +
    scale_fill_manual(values = c("Unique" = "steelblue", "Shared" = "orange")) +
    labs(x = "", y = "Number of genes", title = "Region‑specific prognostic genes") +
    theme_minimal()
  ggsave(file.path(FIG_DIR, "overlap_barplot.png"), p_bar, width = 6, height = 5, dpi = 300)
} else {
  cat("  Not enough regions with significant genes for overlap analysis.\n")
}

# B7. Pairwise HR scatter plots (for regions with common genes)
cat("\n[B7] Pairwise hazard ratio scatter plots...\n")
if (length(region_results) >= 2) {
  reg_names <- names(region_results)
  for (i in 1:(length(reg_names)-1)) {
    for (j in (i+1):length(reg_names)) {
      reg1 <- reg_names[i]; reg2 <- reg_names[j]
      res1 <- region_results[[reg1]]$res; res2 <- region_results[[reg2]]$res
      common_genes <- intersect(res1$Gene, res2$Gene)
      if (length(common_genes) > 0) {
        merged <- merge(res1[, c("Gene", "HR", "FDR")], res2[, c("Gene", "HR", "FDR")],
                        by = "Gene", suffixes = c(paste0("_", gsub(" ", "_", reg1)),
                                                  paste0("_", gsub(" ", "_", reg2))))
        hr_col1 <- paste0("HR_", gsub(" ", "_", reg1))
        hr_col2 <- paste0("HR_", gsub(" ", "_", reg2))
        fdr_col1 <- paste0("FDR_", gsub(" ", "_", reg1))
        fdr_col2 <- paste0("FDR_", gsub(" ", "_", reg2))
        merged$signif_both <- ifelse(merged[[fdr_col1]] < FDR_THRESHOLD & merged[[fdr_col2]] < FDR_THRESHOLD,
                                     "Both sig", "Other")
        p_scatter <- ggplot(merged, aes(x = log2(.data[[hr_col1]]), y = log2(.data[[hr_col2]]), color = signif_both)) +
          geom_point(alpha = 0.6, size = 1.5) + geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
          scale_color_manual(values = c("Both sig" = "red", "Other" = "gray")) +
          labs(x = paste("log2 HR (", reg1, ")"), y = paste("log2 HR (", reg2, ")"), title = "Comparison of effect sizes") +
          theme_classic()
        safe_name <- paste0("HR_scatter_", gsub(" ", "_", reg1), "_vs_", gsub(" ", "_", reg2), ".png")
        ggsave(file.path(FIG_DIR, safe_name), p_scatter, width = 6, height = 5, dpi = 300)
      }
    }
  }
}

# B8. Hallmark pathway enrichment (using msigdbr if available)
cat("\n[B8] Pathway enrichment analysis (Hallmark gene sets)...\n")
if (require("msigdbr", quietly = TRUE)) {
  hallmark_df <- msigdbr(species = "Homo sapiens", category = "H")
  hallmark_list <- split(hallmark_df$gene_symbol, hallmark_df$gs_name)
  cat("  Using official Hallmark gene sets from msigdbr (", length(hallmark_list), " pathways).\n")
} else {
  # Fallback manually defined sets (expanded from original)
  cat("  msigdbr not available. Using built‑in Hallmark‑like sets (less comprehensive).\n")
  hallmark_list <- list(
    "TNFA_NFKB" = c("CCL2","CXCL1","CXCL2","ICAM1","NFKBIA","PTGS2","VCAM1","JUN","FOS","MYC"),
    "HYPOXIA" = c("VEGFA","LDHA","PDK1","SLC2A1","CA9","ENO2","BNIP3","HIF1A","EPO"),
    "EMT" = c("CDH2","VIM","FN1","MMP2","MMP9","SNAI1","SNAI2","TWIST1","ZEB1","CDH1"),
    "APOPTOSIS" = c("BAX","BCL2","CASP3","CASP8","CASP9","FAS","FASLG","BAD","BID"),
    "INFLAMMATORY" = c("IL1B","IL6","TNF","CXCL8","PTGS2","IL10","IL4"),
    "ANGIOGENESIS" = c("VEGFA","VEGFB","VEGFC","KDR","FLT1","PDGFB","ANGPT2","TEK"),
    "MYC_V1" = c("MYC","CCND1","CDK4","ID2","MAX","MXI1","NCL","NPM1"),
    "MYC_V2" = c("CCNA2","CCNB1","CCND2","CDK1","E2F1","MCM2","MCM3","PCNA"),
    "P53" = c("CDKN1A","MDM2","SERPINE1","PMAIP1","BBC3","GADD45A"),
    "G2M" = c("CCNB1","CCNB2","CDK1","PLK1","BUB1","BUB1B"),
    "MITOTIC" = c("AURKA","AURKB","BIRC5","CENPE","KIF11","KIF2C"),
    "MTORC1" = c("RPS6","RPS6KB1","EIF4EBP1","SREBF1","FASN","ACLY")
  )
}

# Background size = number of genes in the expression matrix (after filtering)
bg_size <- ncol(train_expr)
cat("  Background gene set size for hypergeometric test:", bg_size, "\n")

for (reg in names(region_results)) {
  sig_genes <- region_results[[reg]]$sig_genes
  if (length(sig_genes) < 5) next
  enrichment <- data.frame(Pathway = names(hallmark_list), Overlap = NA, p = NA)
  for (i in 1:length(hallmark_list)) {
    overlap <- length(intersect(sig_genes, hallmark_list[[i]]))
    pathway_size <- length(hallmark_list[[i]])
    # Contingency table: overlap, not_overlap_in_pathway, not_overlap_in_background, rest
    cont <- matrix(c(overlap,
                     pathway_size - overlap,
                     length(sig_genes) - overlap,
                     bg_size - length(sig_genes) - (pathway_size - overlap)),
                   nrow = 2)
    ft <- fisher.test(cont, alternative = "greater")
    enrichment$Overlap[i] <- overlap
    enrichment$p[i] <- ft$p.value
  }
  enrichment$FDR <- p.adjust(enrichment$p, method = "BH")
  enrichment <- enrichment[order(enrichment$p), ]
  write.csv(enrichment, file.path(OUTPUT_DIR, paste0(reg, "_enrichment.csv")), row.names = FALSE)

  sig_path <- enrichment[enrichment$FDR < 0.1 & enrichment$Overlap > 0, ]
  if (nrow(sig_path) > 0) {
    sig_path <- head(sig_path, 10)
    sig_path$logP <- -log10(sig_path$p)
    p_path <- ggplot(sig_path, aes(x = logP, y = reorder(Pathway, logP))) +
      geom_bar(stat = "identity", fill = "darkgreen") +
      labs(x = "-log10(p-value)", y = "", title = paste(reg, "- Hallmark enrichment")) +
      theme_minimal()
    ggsave(file.path(FIG_DIR, paste0(gsub(" ", "_", reg), "_enrichment.png")), p_path, width = 7, height = 4, dpi = 300)
  }
}

# B9. Heatmap of common prognostic genes across regions
cat("\n[B9] Heatmap of common prognostic genes...\n")
common_genes_all <- Reduce(intersect, lapply(region_results, function(x) x$sig_genes))
if (length(common_genes_all) > 0) {
  common_ids <- sapply(common_genes_all, function(sym) {
    if (sym %in% names(symbol_to_id)) symbol_to_id[[sym]] else NA
  })
  common_ids <- common_ids[!is.na(common_ids)]
  if (length(common_ids) > 0) {
    top_common_ids <- common_ids[1:min(50, length(common_ids))]
    top_common_syms <- common_genes_all[match(top_common_ids, common_ids)]
    expr_common <- train_expr[, intersect(colnames(train_expr), top_common_ids), drop = FALSE]
    expr_common <- expr_common[rownames(clin_sub), ]
    expr_scaled <- t(scale(expr_common))
    rownames(expr_scaled) <- top_common_syms
    annotation <- data.frame(Subregion = clin_sub$Subregion)
    rownames(annotation) <- rownames(clin_sub)
    pheatmap(expr_scaled, annotation_col = annotation,
             main = paste("Common prognostic genes (n =", length(top_common_syms), ")"),
             fontsize_row = 8, show_colnames = FALSE,
             filename = file.path(FIG_DIR, "common_genes_heatmap.png"), width = 10, height = 10)
  } else {
    cat("  No common genes could be mapped to expression IDs.\n")
  }
} else {
  cat("  No common prognostic genes across regions.\n")
}

# B10. Collect all p‑values for summary table
cat("\n[B10] Collecting p‑values for summary table...\n")
pvalue_table <- data.frame(Analysis = character(), Comparison = character(),
                           P_Value = character(), Significance = character(), stringsAsFactors = FALSE)

for (reg in names(region_results)) {
  res <- region_results[[reg]]$res
  sig_count <- sum(res$FDR < FDR_THRESHOLD, na.rm = TRUE)
  pvalue_table <- rbind(pvalue_table, data.frame(
    Analysis = "Region‑specific Cox",
    Comparison = paste(reg, "- significant genes (FDR<0.05)"),
    P_Value = as.character(sig_count), Significance = ""
  ))
  if (sig_count > 0) {
    top <- res[res$FDR < FDR_THRESHOLD, ][1, ]
    pvalue_table <- rbind(pvalue_table, data.frame(
      Analysis = "Region‑specific Cox",
      Comparison = paste(reg, "- top gene", top$Gene),
      P_Value = format(top$p, scientific = TRUE, digits = 3),
      Significance = ifelse(top$p < 0.001, "***", ifelse(top$p < 0.01, "**", ifelse(top$p < 0.05, "*", "ns")))
    ))
  }
}
# Add subregion log‑rank p‑value from Part A (if available)
if (exists("cox_sub") && !is.null(cox_sub)) {
  logrank_p <- summary(cox_sub)$logtest["pvalue"]
  pvalue_table <- rbind(pvalue_table, data.frame(
    Analysis = "Country analysis", Comparison = "Subregion log‑rank (overall)",
    P_Value = format(logrank_p, scientific = TRUE, digits = 3),
    Significance = ifelse(logrank_p < 0.001, "***", ifelse(logrank_p < 0.01, "**", ifelse(logrank_p < 0.05, "*", "ns")))
  ))
}
write.csv(pvalue_table, file.path(OUTPUT_DIR, "all_pvalues_summary.csv"), row.names = FALSE)

# B11. Pairwise summary table (overlap, HR correlation) – optional
if (length(region_results) >= 2) {
  reg_names <- names(region_results)
  pairwise_summary <- data.frame()
  for (i in 1:(length(reg_names)-1)) {
    for (j in (i+1):length(reg_names)) {
      reg1 <- reg_names[i]; reg2 <- reg_names[j]
      res1 <- region_results[[reg1]]$res; res2 <- region_results[[reg2]]$res
      common_genes <- intersect(res1$Gene, res2$Gene)
      if (length(common_genes) > 0) {
        merged <- merge(res1[, c("Gene", "HR")], res2[, c("Gene", "HR")], by = "Gene")
        merged$logHR1 <- log2(merged$HR.x); merged$logHR2 <- log2(merged$HR.y)
        cor_test <- cor.test(merged$logHR1, merged$logHR2, method = "spearman", use = "complete.obs")
        sig1 <- region_results[[reg1]]$sig_genes; sig2 <- region_results[[reg2]]$sig_genes
        sig_both <- length(intersect(sig1, sig2))
        jaccard <- ifelse(length(union(sig1, sig2)) > 0,
                          length(intersect(sig1, sig2)) / length(union(sig1, sig2)), 0)
        pairwise_summary <- rbind(pairwise_summary, data.frame(
          Region1 = reg1, Region2 = reg2, Common_genes = length(common_genes),
          Correlation_HR_log2 = cor_test$estimate, Correlation_p = cor_test$p.value,
          Significant_in_both = sig_both, Jaccard_index = jaccard
        ))
      }
    }
  }
  if (nrow(pairwise_summary) > 0) {
    write.csv(pairwise_summary, file.path(OUTPUT_DIR, "pairwise_summary_table.csv"), row.names = FALSE)
  }
}

# -----------------------------------------------------------------------------
# FINAL SUMMARY
# -----------------------------------------------------------------------------
cat("\n", rep("=", 80), "\n")
cat("STEP 3 COMPLETED SUCCESSFULLY\n")
cat(rep("=", 80), "\n\n")
cat("Output directory:", OUTPUT_DIR, "\n")
cat("  - Country maps, HR results, KM plots\n")
cat("  - Region‑specific Cox results (CSV), volcano plots, overlap analysis\n")
cat("  - Pathway enrichment, heatmap of common genes\n")
cat("  - All p‑values summary: all_pvalues_summary.csv\n")
cat("Figures saved in:", FIG_DIR, "\n")
cat("End time:", date(), "\n")

sink()
