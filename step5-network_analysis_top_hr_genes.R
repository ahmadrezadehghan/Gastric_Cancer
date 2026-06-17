#!/usr/bin/env Rscript
# =============================================================================
# NETWORK ANALYSIS OF TOP PROGNOSTIC GENES (POSITIVE vs NEGATIVE HR) – ENHANCED
# =============================================================================
# This script:
#   1. Reads univariate Cox results for a chosen region.
#   2. Filters by FDR (optional) and selects top 150 genes with HR > 1 and top 150 with HR < 1.
#   3. Queries STRING database for protein-protein interactions (confidence ≥ 0.4).
#   4. Builds an igraph network, computes degree, and visualizes with ggplot2.
#   5. Saves node/edge tables, unmapped genes list, and publication‑ready figure.
#
# Usage: Rscript network_analysis_top_hr_genes_enhanced.R
# =============================================================================

# -----------------------------------------------------------------------------
# 0. Configuration – EDIT THESE PATHS
# -----------------------------------------------------------------------------
BASE_DIR <- "E:/GastricCancer-2026"
REGION <- "Eastern Asia"   # Change to your region (must match Cox result file)
COX_RESULT_FILE <- file.path(BASE_DIR, "Country_Region_Analysis", paste0(REGION, "_cox_results.csv"))
OUTPUT_DIR <- file.path(BASE_DIR, "Network_Analysis")
STRING_API_KEY <- NULL   # Optional: get from https://string-db.org/help/api/

# Network parameters
N_TOP_POSITIVE <- 150
N_TOP_NEGATIVE <- 150
STRING_CONFIDENCE <- 0.4      # medium confidence (0.4 = 400 on STRING's 0-1000 scale)
FDR_THRESHOLD <- 0.05          # set to NULL to skip FDR filtering
FORCE_REDOWNLOAD <- FALSE

# -----------------------------------------------------------------------------
# 1. Load required packages (install if missing)
# -----------------------------------------------------------------------------
required_pkgs <- c("dplyr", "ggplot2", "igraph", "ggrepel", "readr", "jsonlite", "httr")
if (!require("STRINGdb", quietly = TRUE)) {
  if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")
  BiocManager::install("STRINGdb")
}
library(STRINGdb)
for (pkg in required_pkgs) {
  if (!require(pkg, character.only = TRUE)) {
    install.packages(pkg, dependencies = TRUE)
    library(pkg, character.only = TRUE)
  }
}

# Create output directory
if (!dir.exists(OUTPUT_DIR)) dir.create(OUTPUT_DIR, recursive = TRUE)

# Start log
log_file <- file.path(OUTPUT_DIR, "network_analysis_log.txt")
sink(log_file, append = FALSE, split = TRUE)

cat(rep("=", 80), "\n")
cat("NETWORK ANALYSIS OF TOP PROGNOSTIC GENES (HR > 1 and HR < 1) – ENHANCED\n")
cat("Region:", REGION, "\n")
cat("Start time:", as.character(Sys.time()), "\n")
cat(rep("=", 80), "\n\n")

# -----------------------------------------------------------------------------
# 2. Helper: check internet connectivity
# -----------------------------------------------------------------------------
has_internet <- function() {
  tryCatch({
    GET("https://string-db.org", timeout(5))$status_code == 200
  }, error = function(e) FALSE)
}

# -----------------------------------------------------------------------------
# 3. Load Cox results and select top genes by HR direction
# -----------------------------------------------------------------------------
cat("[1] Loading Cox results from", COX_RESULT_FILE, "\n")
if (!file.exists(COX_RESULT_FILE)) {
  stop("Cox result file not found. Run Step 3 first for region: ", REGION)
}
cox_df <- read.csv(COX_RESULT_FILE, stringsAsFactors = FALSE)

# Keep only genes with valid HR, p, FDR, and gene symbol
cox_df <- cox_df[!is.na(cox_df$HR) & !is.na(cox_df$p) & !is.na(cox_df$FDR) & cox_df$Gene != "", ]
cat("  Total genes with valid HR and FDR:", nrow(cox_df), "\n")

# Optional FDR filter
if (!is.null(FDR_THRESHOLD)) {
  cox_df <- cox_df[cox_df$FDR < FDR_THRESHOLD, ]
  cat("  Genes with FDR <", FDR_THRESHOLD, ":", nrow(cox_df), "\n")
}

if (nrow(cox_df) == 0) stop("No genes passed FDR filtering. Relax threshold or check Cox results.")

# Separate positive and negative HR
pos_hr <- cox_df[cox_df$HR > 1, ]
neg_hr <- cox_df[cox_df$HR < 1, ]

# Sort by HR magnitude and pick top N
pos_hr_sorted <- pos_hr[order(pos_hr$HR, decreasing = TRUE), ]
neg_hr_sorted <- neg_hr[order(neg_hr$HR, decreasing = FALSE), ]

top_pos <- head(pos_hr_sorted$Gene, N_TOP_POSITIVE)
top_neg <- head(neg_hr_sorted$Gene, N_TOP_NEGATIVE)
top_genes <- unique(c(top_pos, top_neg))
cat("  Selected positive HR genes (HR > 1):", length(top_pos), "\n")
cat("  Selected negative HR genes (HR < 1):", length(top_neg), "\n")
cat("  Total unique genes for network:", length(top_genes), "\n")

# Create mapping of gene → HR direction (for node colour)
gene_info <- data.frame(
  Gene = c(top_pos, top_neg),
  HR_direction = c(rep("Positive (HR>1)", length(top_pos)),
                   rep("Negative (HR<1)", length(top_neg))),
  stringsAsFactors = FALSE
)
gene_info <- gene_info[!duplicated(gene_info$Gene), ]

# -----------------------------------------------------------------------------
# 4. Check internet and query STRING
# -----------------------------------------------------------------------------
if (!has_internet()) {
  warning("No internet connection. Cannot query STRING database. Exiting gracefully.")
  # Create a placeholder plot
  placeholder <- ggplot() +
    annotate("text", x = 0.5, y = 0.5, label = "Network not available\n(no internet connection)",
             size = 6, hjust = 0.5) +
    theme_void() +
    labs(title = paste("Protein‑Protein Interaction Network –", REGION))
  ggsave(file.path(OUTPUT_DIR, paste0(gsub(" ", "_", REGION), "_network_no_internet.png")),
         placeholder, width = 8, height = 6, dpi = 150)
  cat("\nNo internet – network generation skipped.\n")
  sink()
  quit(save = "no", status = 0)
}

cat("\n[2] Querying STRING database for interactions...\n")
# Initialize STRINGdb object (human – species 9606)
string_db <- STRINGdb$new(
  version = "11.5",
  species = 9606,
  score_threshold = STRING_CONFIDENCE * 1000,
  input_directory = OUTPUT_DIR,
  api_key = STRING_API_KEY
)

# Map gene symbols to STRING identifiers
cat("  Mapping gene symbols to STRING IDs...\n")
mapping_result <- string_db$map(data.frame(gene = top_genes), "gene", removeUnmappedRows = FALSE)
unmapped <- mapping_result[is.na(mapping_result$STRING_id), "gene"]
mapped <- mapping_result[!is.na(mapping_result$STRING_id), ]
cat("  Successfully mapped:", nrow(mapped), "genes\n")
if (length(unmapped) > 0) {
  cat("  Unmapped genes:", length(unmapped), "\n")
  writeLines(unmapped, file.path(OUTPUT_DIR, paste0(gsub(" ", "_", REGION), "_unmapped_genes.txt")))
}

if (nrow(mapped) == 0) {
  warning("No genes could be mapped to STRING identifiers. Exiting gracefully.")
  placeholder <- ggplot() +
    annotate("text", x = 0.5, y = 0.5, label = "Network not available\n(no genes mapped to STRING)",
             size = 6, hjust = 0.5) +
    theme_void() +
    labs(title = paste("Protein‑Protein Interaction Network –", REGION))
  ggsave(file.path(OUTPUT_DIR, paste0(gsub(" ", "_", REGION), "_network_no_mapping.png")),
         placeholder, width = 8, height = 6, dpi = 150)
  sink()
  quit(save = "no", status = 0)
}

# Get interactions
cat("  Fetching protein-protein interactions (confidence >= ", STRING_CONFIDENCE, ")...\n")
string_interactions <- tryCatch(
  string_db$get_interactions(mapped$STRING_id),
  error = function(e) NULL
)

if (is.null(string_interactions) || nrow(string_interactions) == 0) {
  warning("No interactions found for the selected genes. Try lowering STRING_CONFIDENCE or increasing gene set.")
  placeholder <- ggplot() +
    annotate("text", x = 0.5, y = 0.5,
             label = paste("Network not available\n(no interactions at confidence ≥", STRING_CONFIDENCE, ")"),
             size = 5, hjust = 0.5) +
    theme_void() +
    labs(title = paste("Protein‑Protein Interaction Network –", REGION))
  ggsave(file.path(OUTPUT_DIR, paste0(gsub(" ", "_", REGION), "_network_no_interactions.png")),
         placeholder, width = 8, height = 6, dpi = 150)
  cat("\nNo interactions – network generation skipped.\n")
  sink()
  quit(save = "no", status = 0)
}
cat("  Retrieved", nrow(string_interactions), "unique interactions\n")

# -----------------------------------------------------------------------------
# 5. Build igraph network and compute node attributes
# -----------------------------------------------------------------------------
cat("\n[3] Building network and computing centrality...\n")
g <- graph_from_data_frame(
  d = string_interactions[, c("from", "to")],
  directed = FALSE,
  vertices = data.frame(name = unique(c(string_interactions$from, string_interactions$to)))
)

# Add gene symbol as vertex attribute
id_to_gene <- setNames(mapped$gene, mapped$STRING_id)
V(g)$gene_symbol <- ifelse(V(g)$name %in% names(id_to_gene), id_to_gene[V(g)$name], V(g)$name)

# Add HR direction attribute
V(g)$HR_direction <- "Other (not in top set)"
for (v in V(g)) {
  gene <- V(g)$gene_symbol[v]
  if (gene %in% gene_info$Gene) {
    V(g)$HR_direction[v] <- gene_info$HR_direction[gene_info$Gene == gene]
  }
}

# Compute degree
V(g)$degree <- degree(g)

# Remove isolated nodes (degree 0) from final graph
g_clean <- delete.vertices(g, which(degree(g) == 0))
cat("  Network after removing isolates:", vcount(g_clean), "nodes,", ecount(g_clean), "edges\n")

if (vcount(g_clean) == 0) {
  warning("All nodes isolated after removing degree-0 vertices. Cannot plot network.")
  sink()
  quit(save = "no", status = 0)
}

# -----------------------------------------------------------------------------
# 6. Layout and prepare data for ggplot
# -----------------------------------------------------------------------------
cat("\n[4] Computing layout (this may take a few seconds)...\n")
set.seed(42)
layout <- layout_with_fr(g_clean, niter = 1000)

edge_df <- as_data_frame(g_clean, what = "edges")
edge_df$from_x <- layout[edge_df$from, 1]
edge_df$from_y <- layout[edge_df$from, 2]
edge_df$to_x <- layout[edge_df$to, 1]
edge_df$to_y <- layout[edge_df$to, 2]

node_df <- data.frame(
  name = V(g_clean)$name,
  gene = V(g_clean)$gene_symbol,
  HR_direction = V(g_clean)$HR_direction,
  degree = V(g_clean)$degree,
  x = layout[, 1],
  y = layout[, 2],
  stringsAsFactors = FALSE
)

# Colour mapping
color_map <- c(
  "Positive (HR>1)" = "#D95F02",
  "Negative (HR<1)" = "#1B9E77",
  "Other (not in top set)" = "#7570B3"
)

# -----------------------------------------------------------------------------
# 7. Plot network with ggplot2 (including degree legend)
# -----------------------------------------------------------------------------
cat("\n[5] Generating publication‑ready network plot...\n")
p_network <- ggplot() +
  geom_segment(data = edge_df,
               aes(x = from_x, y = from_y, xend = to_x, yend = to_y),
               color = "gray80", linewidth = 0.3, alpha = 0.6) +
  geom_point(data = node_df,
             aes(x = x, y = y, fill = HR_direction, size = degree),
             shape = 21, color = "black", stroke = 0.3, alpha = 0.9) +
  geom_text_repel(data = node_df %>% arrange(desc(degree)) %>% head(20),
                  aes(x = x, y = y, label = gene),
                  size = 3, max.overlaps = 20, box.padding = 0.5, point.padding = 0.2) +
  scale_fill_manual(values = color_map, name = "Hazard Ratio direction") +
  scale_size_continuous(range = c(1, 8), name = "Degree (number of interactions)") +
  guides(size = guide_legend(title = "Degree", override.aes = list(shape = 21, fill = "gray"))) +
  theme_void() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
    plot.subtitle = element_text(hjust = 0.5, size = 10, color = "gray40"),
    legend.position = "right",
    legend.box = "vertical",
    legend.title = element_text(face = "bold")
  ) +
  labs(
    title = paste("Protein‑Protein Interaction Network of Top Prognostic Genes –", REGION),
    subtitle = paste0(
      "Positive HR (", length(top_pos), " genes) | Negative HR (", length(top_neg),
      " genes) | Nodes: ", nrow(node_df), " | Edges: ", nrow(edge_df)
    )
  )

# Save
png_file <- file.path(OUTPUT_DIR, paste0(gsub(" ", "_", REGION), "_network.png"))
pdf_file <- file.path(OUTPUT_DIR, paste0(gsub(" ", "_", REGION), "_network.pdf"))
ggsave(png_file, p_network, width = 14, height = 12, dpi = 300, bg = "white")
ggsave(pdf_file, p_network, width = 14, height = 12, dpi = 300, bg = "white")
cat("  Network plot saved:\n  -", png_file, "\n  -", pdf_file, "\n")

# -----------------------------------------------------------------------------
# 8. Export node/edge tables and unmapped genes
# -----------------------------------------------------------------------------
cat("\n[6] Exporting node and edge tables (CSV)...\n")
write.csv(node_df, file.path(OUTPUT_DIR, paste0(gsub(" ", "_", REGION), "_network_nodes.csv")),
          row.names = FALSE)
write.csv(edge_df, file.path(OUTPUT_DIR, paste0(gsub(" ", "_", REGION), "_network_edges.csv")),
          row.names = FALSE)
cat("  Node table:", nrow(node_df), "rows\n")
cat("  Edge table:", nrow(edge_df), "rows\n")

# -----------------------------------------------------------------------------
# 9. Network summary statistics
# -----------------------------------------------------------------------------
cat("\n[7] Network summary statistics\n")
cat("  - Nodes (proteins):", vcount(g_clean), "\n")
cat("  - Edges (interactions):", ecount(g_clean), "\n")
cat("  - Density:", round(edge_density(g_clean), 4), "\n")
cat("  - Average degree:", round(mean(degree(g_clean)), 2), "\n")
cat("  - Connected components:", components(g_clean)$no, "\n")

summary_stats <- list(
  region = REGION,
  fdr_threshold = FDR_THRESHOLD,
  n_top_positive = length(top_pos),
  n_top_negative = length(top_neg),
  n_unique_genes = length(top_genes),
  n_string_mapped = nrow(mapped),
  n_unmapped = length(unmapped),
  n_network_nodes = vcount(g_clean),
  n_network_edges = ecount(g_clean),
  network_density = edge_density(g_clean),
  average_degree = mean(degree(g_clean)),
  components = components(g_clean)$no,
  string_confidence_threshold = STRING_CONFIDENCE,
  date = as.character(Sys.time())
)
write_json(summary_stats, file.path(OUTPUT_DIR, paste0(gsub(" ", "_", REGION), "_network_summary.json")),
           pretty = TRUE, auto_unbox = TRUE)

cat("\n", rep("=", 80), "\n")
cat("NETWORK ANALYSIS COMPLETED SUCCESSFULLY\n")
cat("Output directory:", OUTPUT_DIR, "\n")
cat("End time:", as.character(Sys.time()), "\n")
sink()
