#!/usr/bin/env Rscript
# =============================================================================
# COMPLETE VISUALIZATIONS: CLINICAL OVERVIEW + SURVIVAL MODEL PERFORMANCE
# =============================================================================
# This script produces publication‑ready figures:
#   1. Clinical data overview (expression dimensions, missingness, age, gender,
#      AJCC stage, T/N stage, KM curve, country, organisation)
#   2. Survival model performance (KM by risk group, risk score distribution,
#      time‑dependent AUC, model comparison C‑index if available)
#
# All figures are saved as high‑resolution PNG/PDF in a dedicated folder.
# =============================================================================

# -----------------------------------------------------------------------------
# 0. User configuration (EDIT THESE PATHS AND PARAMETERS)
# -----------------------------------------------------------------------------
BASE_DIR <- "E:/GastricCancer-2026"
PLATFORM <- "rnaseq"   # either "rnaseq" or "microarray" – must match Step 1.1 output
CLINICAL_FILE <- file.path(BASE_DIR, "Clinical_Data", "demo.xlsx")   # same as Step 2
EXPR_TRAIN_FILE <- file.path(BASE_DIR, "Prepared_Data", PLATFORM, "train_data.csv")
EXPR_TEST_FILE  <- file.path(BASE_DIR, "Prepared_Data", PLATFORM, "test_data.csv")
RISK_SCORES_FILE <- file.path(BASE_DIR, "Survival_ElasticNet_DCA", "test_risk_scores.csv")
SUMMARY_JSON <- file.path(BASE_DIR, "Survival_ElasticNet_DCA", "summary.json")
OUTPUT_FIG_DIR <- file.path(BASE_DIR, "Manuscript_Figures")

# Create output directory
if (!dir.exists(OUTPUT_FIG_DIR)) dir.create(OUTPUT_FIG_DIR, recursive = TRUE)

# -----------------------------------------------------------------------------
# 1. Load required packages (install if missing)
# -----------------------------------------------------------------------------
required_pkgs <- c("readxl", "ggplot2", "survival", "survminer", "patchwork",
                   "viridis", "tidyverse", "jsonlite", "gridExtra")
for (pkg in required_pkgs) {
  if (!require(pkg, character.only = TRUE)) {
    install.packages(pkg, dependencies = TRUE)
    library(pkg, character.only = TRUE)
  }
}

# Set global theme for publication
theme_set(theme_minimal(base_size = 12) +
            theme(
              plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
              plot.subtitle = element_text(hjust = 0.5, size = 10, color = "gray40"),
              legend.position = "bottom",
              panel.grid.minor = element_blank(),
              panel.border = element_rect(fill = NA, color = "gray80", linewidth = 0.5),
              axis.text.x = element_text(angle = 45, hjust = 1)
            ))

# -----------------------------------------------------------------------------
# PART 1: CLINICAL DATA OVERVIEW FIGURE (Panels A-J)
# -----------------------------------------------------------------------------
cat("\n", rep("=", 80), "\n")
cat("PART 1: CLINICAL DATA OVERVIEW\n")
cat(rep("=", 80), "\n")

# 1.1 Load clinical data (same as Step 2)
cat("Loading clinical data from", CLINICAL_FILE, "...\n")
if (!file.exists(CLINICAL_FILE)) stop("Clinical file not found: ", CLINICAL_FILE)
demo_raw <- read_excel(CLINICAL_FILE, col_names = TRUE)
features <- demo_raw[[1]]
demo_mat <- as.matrix(demo_raw[, -1])
rownames(demo_mat) <- features
colnames(demo_mat) <- colnames(demo_raw)[-1]

# Helper to parse decimal commas (same as Step 2)
parse_decimal <- function(x) {
  if (length(x) == 0 || is.na(x) || x == "") return(NA_real_)
  x <- as.character(x)
  x <- gsub(",", ".", x)
  x <- trimws(x)
  as.numeric(x)
}

# Extract clinical variables (same columns as Step 2)
os_row    <- grep("overall_survival", rownames(demo_mat), ignore.case = TRUE)[1]
status_row <- grep("survival status", rownames(demo_mat), ignore.case = TRUE)[1]
stage_row <- grep("^stage$", rownames(demo_mat), ignore.case = TRUE)[1]
nstage_row <- grep("n_stage", rownames(demo_mat), ignore.case = TRUE)[1]
age_row   <- grep("^age$", rownames(demo_mat), ignore.case = TRUE)[1]
gender_row <- grep("^gender$", rownames(demo_mat), ignore.case = TRUE)[1]
country_row <- grep("country", rownames(demo_mat), ignore.case = TRUE)[1]
org_row    <- grep("organization", rownames(demo_mat), ignore.case = TRUE)[1]

sample_ids <- colnames(demo_mat)
clinical_list <- list()
for (sid in sample_ids) {
  clinical_list[[sid]] <- data.frame(
    sample_id = sid,
    Age = parse_decimal(demo_mat[age_row, sid]),
    Gender = as.character(demo_mat[gender_row, sid]),
    Stage = as.character(demo_mat[stage_row, sid]),
    N_stage = parse_decimal(demo_mat[nstage_row, sid]),
    OS_time = parse_decimal(demo_mat[os_row, sid]),
    OS_status = ifelse(tolower(demo_mat[status_row, sid]) %in% c("dead", "deceased"), 1, 0),
    Country = as.character(demo_mat[country_row, sid]),
    Organization = as.character(demo_mat[org_row, sid]),
    stringsAsFactors = FALSE
  )
}
clinical_data <- do.call(rbind, clinical_list)
rownames(clinical_data) <- clinical_data$sample_id
clinical_data$sample_id <- NULL
clinical_data <- clinical_data[!is.na(clinical_data$OS_time) & !is.na(clinical_data$OS_status) & clinical_data$OS_time > 0, ]
cat("  Samples in clinical data:", nrow(clinical_data), "\n")

# 1.2 Expression dimensions (train + test)
cat("Loading expression data...\n")
if (!file.exists(EXPR_TRAIN_FILE)) stop("Train expression file not found: ", EXPR_TRAIN_FILE)
if (!file.exists(EXPR_TEST_FILE)) stop("Test expression file not found: ", EXPR_TEST_FILE)
train_expr <- read.csv(EXPR_TRAIN_FILE, row.names = 1, check.names = FALSE)
test_expr <- read.csv(EXPR_TEST_FILE, row.names = 1, check.names = FALSE)
total_samples <- nrow(train_expr) + nrow(test_expr)
total_genes <- ncol(train_expr)
cat("  Expression matrix:", total_genes, "genes ×", total_samples, "samples\n")

# 1.3 Create individual panels
# Panel A: Expression matrix dimensions
panel_a <- ggplot() +
  annotate("text", x = 0.5, y = 0.8, size = 6, fontface = "bold",
           label = paste0("Expression matrix\n", total_genes, " genes × ", total_samples, " samples")) +
  annotate("text", x = 0.5, y = 0.5, size = 4,
           label = paste0("Train: ", nrow(train_expr), " samples\nTest: ", nrow(test_expr), " samples")) +
  theme_void() + ggtitle("A) Expression data") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))

# Panel B: Missingness (bar plot)
miss_df <- data.frame(Variable = names(clinical_data),
                      Missing = sapply(clinical_data, function(x) sum(is.na(x)) / nrow(clinical_data)))
miss_df <- miss_df[miss_df$Missing > 0, ]
if (nrow(miss_df) == 0) miss_df <- data.frame(Variable = "None", Missing = 0)
panel_b <- ggplot(miss_df, aes(x = reorder(Variable, -Missing), y = Missing, fill = Variable)) +
  geom_bar(stat = "identity") + scale_y_continuous(labels = scales::percent) +
  labs(title = "B) Missing proportions", x = "", y = "% missing") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"), legend.position = "none")

# Panel C: Age distribution
panel_c <- ggplot(clinical_data, aes(x = Age)) +
  geom_histogram(fill = "steelblue", color = "black", bins = 30, na.rm = TRUE) +
  labs(title = "C) Age distribution", x = "Age (years)", y = "Count") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))

# Panel D: Gender (only Male/Female)
gender_counts <- table(clinical_data$Gender)
gender_df <- data.frame(Gender = names(gender_counts), Count = as.numeric(gender_counts))
panel_d <- ggplot(gender_df, aes(x = Gender, y = Count, fill = Gender)) +
  geom_bar(stat = "identity") + scale_fill_manual(values = c("Male" = "#2c7bb6", "Female" = "#fdae61")) +
  labs(title = "D) Gender", x = "", y = "Count") + theme(legend.position = "none")

# Panel E: AJCC Stage
stage_counts <- table(clinical_data$Stage)
stage_df <- data.frame(Stage = names(stage_counts), Count = as.numeric(stage_counts))
stage_df <- stage_df[!is.na(stage_df$Stage), ]
panel_e <- ggplot(stage_df, aes(x = Stage, y = Count, fill = Stage)) +
  geom_bar(stat = "identity") + scale_fill_viridis_d() +
  labs(title = "E) AJCC Stage", x = "Stage", y = "Count") + theme(legend.position = "none")

# Panel F: N stage (numeric)
n_df <- clinical_data[!is.na(clinical_data$N_stage), ]
panel_f <- ggplot(n_df, aes(x = factor(N_stage))) +
  geom_bar(fill = "steelblue") +
  labs(title = "F) N stage", x = "N stage", y = "Count") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))

# Panel G: Overall survival Kaplan‑Meier
km_fit <- survfit(Surv(OS_time, OS_status) ~ 1, data = clinical_data)
panel_g <- ggsurvplot(km_fit, data = clinical_data, risk.table = FALSE, conf.int = TRUE,
                      palette = "black", title = "G) Overall survival",
                      xlab = "Time (months)", ylab = "Survival probability", legend = "none")$plot +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))

# Panel H: Country (top 5)
country_tab <- sort(table(clinical_data$Country), decreasing = TRUE)
country_df <- data.frame(Country = names(country_tab)[1:min(5, length(country_tab))],
                         Count = as.numeric(country_tab)[1:min(5, length(country_tab))])
panel_h <- ggplot(country_df, aes(x = reorder(Country, -Count), y = Count, fill = Country)) +
  geom_bar(stat = "identity") + scale_fill_viridis_d() +
  labs(title = "H) Country (top 5)", x = "", y = "Count") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none")

# Panel I: Organization (top 6)
org_tab <- sort(table(clinical_data$Organization), decreasing = TRUE)
org_df <- data.frame(Organization = names(org_tab)[1:min(6, length(org_tab))],
                     Count = as.numeric(org_tab)[1:min(6, length(org_tab))])
panel_i <- ggplot(org_df, aes(x = reorder(Organization, -Count), y = Count, fill = Organization)) +
  geom_bar(stat = "identity") + scale_fill_viridis_d() +
  labs(title = "I) Organization (top 6)", x = "", y = "Count") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none")

# Combine clinical panels (now A-I, 9 panels)
clinical_combined <- (panel_a + panel_b) /
  (panel_c + panel_d + panel_e) /
  (panel_f + panel_g + panel_h) /
  (panel_i + plot_spacer()) +
  plot_annotation(title = "Clinical and Expression Data Overview",
                  theme = theme(plot.title = element_text(size = 18, face = "bold", hjust = 0.5))) +
  plot_layout(heights = c(0.8, 1, 1, 0.8))

# Save clinical figure
ggsave(file.path(OUTPUT_FIG_DIR, "Figure_S1_Clinical_Overview.png"),
       clinical_combined, width = 16, height = 18, dpi = 300, bg = "white")
cat("  Saved: Figure_S1_Clinical_Overview.png\n")

# -----------------------------------------------------------------------------
# PART 2: SURVIVAL MODEL PERFORMANCE FIGURES
# -----------------------------------------------------------------------------
cat("\n", rep("=", 80), "\n")
cat("PART 2: SURVIVAL MODEL PERFORMANCE\n")
cat(rep("=", 80), "\n")

# 2.1 Load risk scores (from Step 2)
if (!file.exists(RISK_SCORES_FILE)) stop("Risk scores file not found: ", RISK_SCORES_FILE)
risk_df <- read.csv(RISK_SCORES_FILE)
risk_df$OS_time <- as.numeric(risk_df$OS_time)
risk_df$OS_status <- as.numeric(risk_df$OS_status)

# Determine risk score column (try common names)
if ("risk_score" %in% names(risk_df)) {
  risk_df$primary_risk <- risk_df$risk_score
} else if ("risk_clingen_optimized" %in% names(risk_df)) {
  risk_df$primary_risk <- risk_df$risk_clingen_optimized
} else {
  stop("No risk score column found. Expected 'risk_score' or 'risk_clingen_optimized'.")
}

# Check for clinical‑only risk column (optional)
has_clin_only <- "risk_clin_only" %in% names(risk_df)

median_risk <- median(risk_df$primary_risk, na.rm = TRUE)
risk_df$risk_group <- factor(ifelse(risk_df$primary_risk > median_risk, "High Risk", "Low Risk"),
                             levels = c("Low Risk", "High Risk"))
cat("  Test samples:", nrow(risk_df), " events:", sum(risk_df$OS_status), "\n")
cat("  Median risk score:", round(median_risk, 4), "\n")

# 2.2 Kaplan‑Meier curves
fit_km <- survfit(Surv(OS_time, OS_status) ~ risk_group, data = risk_df)
sdiff <- survdiff(Surv(OS_time, OS_status) ~ risk_group, data = risk_df)
p_val <- 1 - pchisq(sdiff$chisq, length(sdiff$n) - 1)
p_text <- ifelse(p_val < 0.001, "p < 0.001", paste("p =", round(p_val, 3)))

km_plot <- ggsurvplot(fit_km, data = risk_df, risk.table = TRUE, pval = FALSE,
                      conf.int = TRUE, palette = c("#2E86AB", "#C73E1D"),
                      xlab = "Time (months)", ylab = "Overall Survival Probability",
                      title = "Kaplan-Meier Survival Curves",
                      legend.title = "Risk Group", legend.labs = c("Low Risk", "High Risk"),
                      risk.table.height = 0.25, ggtheme = theme_minimal())
km_plot_obj <- km_plot$plot +
  annotate("text", x = max(risk_df$OS_time, na.rm = TRUE) * 0.7, y = 0.9,
           label = p_text, size = 4, hjust = 0) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))

# 2.3 Risk score distribution
risk_colors <- c("Low Risk" = "#2E86AB", "High Risk" = "#C73E1D")
risk_dist_plot <- ggplot(risk_df, aes(x = primary_risk, fill = risk_group, color = risk_group)) +
  geom_density(alpha = 0.5, linewidth = 0.8) +
  scale_fill_manual(values = risk_colors) + scale_color_manual(values = risk_colors) +
  labs(title = "Risk Score Distribution",
       subtitle = paste("High Risk (n =", sum(risk_df$risk_group == "High Risk"),
                        "), Low Risk (n =", sum(risk_df$risk_group == "Low Risk"), ")"),
       x = "Risk Score", y = "Density", fill = "Risk Group", color = "Risk Group")

# 2.4 Time‑dependent AUC from Step 2 summary.json
if (!file.exists(SUMMARY_JSON)) stop("Summary JSON not found: ", SUMMARY_JSON)
summary_json <- jsonlite::read_json(SUMMARY_JSON, simplifyVector = TRUE)

# Extract time‑dependent AUC list (each element: time, auc, ci_lower, ci_upper)
auc_list <- summary_json$time_dependent_auc
if (is.null(auc_list)) {
  stop("No time_dependent_auc in summary.json. Check Step 2 output.")
}
auc_df <- do.call(rbind, lapply(auc_list, function(x) {
  data.frame(Time = as.numeric(x$time), AUC = as.numeric(x$auc),
             CI_lower = as.numeric(x$ci_lower), CI_upper = as.numeric(x$ci_upper))
}))
auc_df$Time <- factor(auc_df$Time, levels = sort(unique(auc_df$Time)))

auc_plot <- ggplot(auc_df, aes(x = Time, y = AUC, group = 1)) +
  geom_line(color = "#2E86AB", linewidth = 1.2) +
  geom_point(size = 4, color = "#C73E1D") +
  geom_errorbar(aes(ymin = CI_lower, ymax = CI_upper), width = 0.15, color = "#C73E1D", linewidth = 0.8) +
  geom_text(aes(label = sprintf("%.3f", AUC)), vjust = -1.5, size = 4) +
  ylim(0.5, 0.85) +
  labs(title = "Time-Dependent AUC",
       subtitle = "From training‑derived model on test set (95% CI)",
       x = "Prediction Horizon (months)", y = "AUC")

# 2.5 Model comparison C‑index (if clinical‑only risk column exists)
if (has_clin_only) {
  c_clin <- summary(coxph(Surv(OS_time, OS_status) ~ risk_clin_only, data = risk_df))$concordance[1]
  c_gene <- summary(coxph(Surv(OS_time, OS_status) ~ primary_risk, data = risk_df))$concordance[1]
  cindex_vals <- data.frame(Model = c("Clinical Only", "Clinical + Gene"), Cindex = c(c_clin, c_gene))
  bar_plot <- ggplot(cindex_vals, aes(x = reorder(Model, Cindex), y = Cindex, fill = Model)) +
    geom_bar(stat = "identity", width = 0.7, alpha = 0.85) +
    geom_text(aes(label = sprintf("%.4f", Cindex)), hjust = -0.1, size = 4, fontface = "bold") +
    geom_hline(yintercept = 0.5, linetype = "dashed", color = "gray50") +
    scale_fill_manual(values = c("#2E86AB", "#C73E1D")) + coord_flip(ylim = c(0, 0.8)) +
    labs(title = "Model Comparison (C-index)",
         subtitle = paste("Improvement: +", round(c_gene - c_clin, 4)),
         x = "", y = "C-index") + theme(legend.position = "none")
  cat("  Clinical‑only risk scores found – including comparison bar plot.\n")
} else {
  cat("  'risk_clin_only' not found in risk scores; skipping model comparison bar plot.\n")
  bar_plot <- ggplot() + annotate("text", x = 0.5, y = 0.5, label = "Clinical-only model not available") + theme_void()
}

# 2.6 Combine model performance figures
model_combined <- (auc_plot + km_plot_obj) / (risk_dist_plot + bar_plot) +
  plot_annotation(
    title = "SURVIVAL MODEL PERFORMANCE SUMMARY",
    subtitle = paste0("AUC (3‑year = ", round(auc_df$AUC[auc_df$Time == 36], 3),
                      "); C-index = ", round(ifelse(has_clin_only, c_gene, summary_json$test_c_index), 4),
                      " | Test set: n = ", nrow(risk_df), ", events = ", sum(risk_df$OS_status)),
    theme = theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
                  plot.subtitle = element_text(hjust = 0.5, size = 12))
  )

ggsave(file.path(OUTPUT_FIG_DIR, "Figure_2_Survival_Performance.png"),
       model_combined, width = 14, height = 12, dpi = 300, bg = "white")
cat("  Saved: Figure_2_Survival_Performance.png\n")

# Also save individual panels for use in manuscripts
ggsave(file.path(OUTPUT_FIG_DIR, "KM_curves.png"), km_plot_obj, width = 7, height = 7, dpi = 300)
ggsave(file.path(OUTPUT_FIG_DIR, "Risk_distribution.png"), risk_dist_plot, width = 7, height = 6, dpi = 300)
ggsave(file.path(OUTPUT_FIG_DIR, "Time_dependent_AUC.png"), auc_plot, width = 6, height = 5, dpi = 300)
if (has_clin_only) {
  ggsave(file.path(OUTPUT_FIG_DIR, "Model_comparison_Cindex.png"), bar_plot, width = 8, height = 5, dpi = 300)
}

# 2.7 Export summary tables
risk_summary <- risk_df %>% group_by(risk_group) %>%
  summarise(N = n(), Events = sum(OS_status), Median_OS = median(OS_time, na.rm = TRUE), .groups = "drop")
write.csv(risk_summary, file.path(OUTPUT_FIG_DIR, "risk_group_summary.csv"), row.names = FALSE)

if (has_clin_only) write.csv(cindex_vals, file.path(OUTPUT_FIG_DIR, "model_comparison_cindex.csv"), row.names = FALSE)
write.csv(auc_df, file.path(OUTPUT_FIG_DIR, "time_dependent_auc_summary.csv"), row.names = FALSE)

# -----------------------------------------------------------------------------
# FINAL SUMMARY
# -----------------------------------------------------------------------------
cat("\n", rep("=", 80), "\n")
cat("ALL VISUALIZATIONS COMPLETED SUCCESSFULLY\n")
cat(rep("=", 80), "\n")
cat("Output directory:", OUTPUT_FIG_DIR, "\n")
cat("Files created:\n")
cat("  - Figure_S1_Clinical_Overview.png (clinical data overview)\n")
cat("  - Figure_2_Survival_Performance.png (combined model figures)\n")
cat("  - KM_curves.png, Risk_distribution.png, Time_dependent_AUC.png")
if (has_clin_only) cat(", Model_comparison_Cindex.png")
cat("\n  - risk_group_summary.csv, time_dependent_auc_summary.csv\n")
cat("End time:", as.character(Sys.time()), "\n")
