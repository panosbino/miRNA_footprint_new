#!/usr/bin/env Rscript

# ==============================================================================
# GFP fluorescence QC: merge FACS fluorescence with raw and scran-normalized
# GFP transgene counts, and compare correlations.
#
# Follow-up to build_annotated_scran_dataset.R -- loads its saved
# cell_metadata.rds (which already has GFP_raw_counts, GFP_normalized,
# Dox_Concentration/Group, and UMAP1/UMAP2) rather than re-running scran
# normalization, since that step is expensive and meant to run once.
#
# Ported from GFP_flu_count_cor.R (Seurat-based). No Seurat used here.
#
# FIXES APPLIED vs. the original:
#   1. Barcode/well matching now uses intersect() instead of directly
#      indexing well_map by filtered_barcodes as row names. The original
#          well_map[filtered_barcodes, c("WellID","bc_set")]
#      silently returns a row of NAs for any filtered barcode not present in
#      the well map, which then flowed into merge() with no warning. Now the
#      match rate is checked and a warning is raised if it looks low, same
#      pattern used for the Dox/barcode merge elsewhere in this project.
#   2. merge() row counts are now logged before/after, so silently dropped
#      GFP-fluorescence or barcode rows (e.g. from a WellID formatting
#      mismatch) are visible instead of just producing a smaller-than-
#      expected result with no explanation.
#   3. Compares BOTH raw GFP counts and scran-normalized GFP against
#      fluorescence (previously only compared fluorescence against Seurat's
#      RC-normalized value, which is a different, non-scran normalization
#      and isn't what "normalized GFP" means in the current pipeline).
#   4. Correlation-label plot annotations are now built from the actual
#      per-group summary table (matched by group name) instead of hardcoded
#      position vectors (vjust = seq(1.5, 4.5, by=1), fontface/color vectors
#      assuming exactly 4 rows in a fixed order). The original would
#      silently mislabel/miscolor rows if a Dox level had zero cells or the
#      group order changed.
#   5. Groups with fewer than MIN_CELLS_FOR_COR cells are flagged and
#      excluded from correlation plots/summaries rather than silently
#      plotted as if reliable.
#   6. Non-positive values are removed (with a logged count) before log10
#      scaling, instead of being silently dropped by ggplot2's log scale
#      with just a generic "removed rows" warning.
#   7. GFP gene name resolved dynamically (reused from the earlier
#      annotation step) instead of hardcoded "eGFP".
# ==============================================================================

## ---- Config -----------------------------------------------------------------
RAW_DATA_DIR <- "/cfs/klemming/projects/snic/naiss2024-6-235/miRNA_footprint/datasets/HEK_SS3/raw"
PROCESSED_DATA_DIR <- "/cfs/klemming/projects/snic/naiss2024-6-235/miRNA_footprint/datasets/HEK_SS3/processed"
PLOT_DIR <- "/cfs/klemming/projects/snic/naiss2024-6-235/miRNA_footprint/analysis/HEK_SS3/qc/plots"
PLATES <- c("101", "103", "105")
SEQUENCING_PROJECT_PREFIX <- "P32156"
MIN_BARCODE_MATCH_FRAC <- 0.5
MIN_CELLS_FOR_COR <- 10          # groups smaller than this are flagged, not trusted
INDUCED_LEVELS <- c("0.01 Dox", "0.1 Dox", "1 Dox")

## ---- Libraries ----------------------------------------------------------------
library(tidyverse)
library(ggplot2)
library(patchwork)
library(ggridges)

## ---- Load annotated cell metadata from build_annotated_scran_dataset.R --
cell_metadata <- readRDS(file.path(PROCESSED_DATA_DIR, "cell_metadata.rds"))

# Resolve GFP gene id the same way as the main pipeline (kept here only for
# labeling plots/messages; the actual raw/normalized GFP values are already
# in cell_metadata).
possible_gfp_names <- c("eGFP", "GFP", "Gfp", "EGFP", "EGfp")
egfp_gene_id <- possible_gfp_names[possible_gfp_names %in% "eGFP"][1]  # matches build script's default naming
if (is.na(egfp_gene_id)) egfp_gene_id <- "GFP"

## ============================================================================
## 1. Load and merge GFP fluorescence (FACS) data per plate
## ============================================================================

process_gfp_fluorescence <- function(plate_id) {
  gfp_file <- file.path(RAW_DATA_DIR, paste0("P", plate_id, "_GFP_fluorescence.tsv"))
  gfp_data <- read.delim(gfp_file, sep = "\t", header = TRUE)

  well_map_file <- file.path(RAW_DATA_DIR,
                              paste0(SEQUENCING_PROJECT_PREFIX, "_", plate_id, ".well_barcodes.txt"))
  well_map <- read.delim(well_map_file, header = TRUE, sep = "\t")
  rownames(well_map) <- well_map$bc_set

  barcode_file <- file.path(RAW_DATA_DIR, paste0("cell_barcodes_filtered_", plate_id, ".txt"))
  filtered_barcodes <- readLines(barcode_file)

  # FIX: intersect first instead of indexing well_map by filtered_barcodes
  # directly -- the latter silently produces NA rows for any barcode not
  # present in well_map, with no warning.
  matched_barcodes <- intersect(filtered_barcodes, rownames(well_map))
  if (length(matched_barcodes) == 0) {
    warning(sprintf("Plate %s: no filtered barcodes matched the well map -- GFP fluorescence will be entirely missing for this plate.", plate_id))
  } else {
    match_frac <- length(matched_barcodes) / length(filtered_barcodes)
    if (match_frac < MIN_BARCODE_MATCH_FRAC) {
      warning(sprintf("Plate %s: only %d / %d filtered barcodes (%.1f%%) matched the well map.",
                       plate_id, length(matched_barcodes), length(filtered_barcodes), 100 * match_frac))
    }
  }
  well_barcode_map <- well_map[matched_barcodes, c("WellID", "bc_set")]

  n_before <- nrow(well_barcode_map)
  merged_data <- merge(gfp_data, well_barcode_map, by = "WellID", all.x = FALSE)
  n_after <- nrow(merged_data)
  message(sprintf(
    "  Plate %s: %d barcoded wells, %d matched a fluorescence WellID (%d dropped -- check WellID formatting if this is unexpectedly high)",
    plate_id, n_before, n_after, n_before - n_after
  ))

  data.frame(
    cell_id = paste0(plate_id, "_", merged_data$bc_set),
    GFP_fluorescence = merged_data$GFP,
    row.names = paste0(plate_id, "_", merged_data$bc_set)
  )
}

message("Processing GFP fluorescence data for each plate...")
all_gfp_data <- do.call(rbind, lapply(PLATES, process_gfp_fluorescence))

cell_metadata$GFP_fluorescence <- NA_real_
matched_cells <- intersect(rownames(all_gfp_data), cell_metadata$cell_id)
rownames(cell_metadata) <- cell_metadata$cell_id
cell_metadata[matched_cells, "GFP_fluorescence"] <- all_gfp_data[matched_cells, "GFP_fluorescence"]

n_missing_flu <- sum(is.na(cell_metadata$GFP_fluorescence))
message(sprintf("%d / %d cells (%.1f%%) have no matched GFP fluorescence value.",
                 n_missing_flu, nrow(cell_metadata), 100 * n_missing_flu / nrow(cell_metadata)))

## ============================================================================
## 2. UMAP / ridge plots colored by GFP fluorescence
## ============================================================================

theme_set(theme_minimal())
flu_99 <- quantile(cell_metadata$GFP_fluorescence, 0.99, na.rm = TRUE)

p_flu_umap <- ggplot(cell_metadata, aes(UMAP1, UMAP2, color = pmin(GFP_fluorescence, flu_99))) +
  geom_point(size = 2, alpha = 0.8) +
  scale_color_distiller(palette = "YlOrRd", direction = 1, name = "GFP\nfluorescence") +
  ggtitle("GFP fluorescence on UMAP")
ggsave(file.path(PLOT_DIR, "umap_gfp_fluorescence.png"), p_flu_umap, width = 10, height = 8, dpi = 300)

p_flu_by_dox <- ggplot(cell_metadata %>% filter(Dox_Concentration != "Unknown"),
                        aes(UMAP1, UMAP2, color = pmin(GFP_fluorescence, flu_99))) +
  geom_point(size = 1.5, alpha = 0.8) +
  facet_wrap(~ Dox_Concentration, nrow = 1) +
  scale_color_distiller(palette = "YlOrRd", direction = 1, name = "GFP\nfluorescence") +
  ggtitle("GFP fluorescence by Dox concentration")
ggsave(file.path(PLOT_DIR, "umap_gfp_fluorescence_by_dox.png"), p_flu_by_dox, width = 24, height = 8, dpi = 300)

flu_ridge <- ggplot(cell_metadata, aes(x = pmin(GFP_fluorescence, flu_99 * 1.1),
                                        y = Dox_Concentration, fill = Dox_Concentration)) +
  geom_density_ridges(alpha = 0.7) +
  scale_fill_brewer(palette = "Dark2") +
  labs(x = "GFP fluorescence", y = "Dox Concentration",
       title = "GFP fluorescence distribution by Dox concentration")
ggsave(file.path(PLOT_DIR, "ridge_gfp_fluorescence_by_dox.png"), flu_ridge, width = 10, height = 6, dpi = 300)

## ============================================================================
## 3. Correlations: fluorescence vs. raw GFP counts and vs. scran-normalized GFP
## ============================================================================

correlation_data <- cell_metadata %>%
  select(cell_id, plate, Dox_Concentration, Dox_Group, GFP_fluorescence, GFP_raw_counts, GFP_normalized)

compute_group_correlations <- function(df, group_var, min_n = MIN_CELLS_FOR_COR) {
  df %>%
    group_by(.data[[group_var]]) %>%
    summarize(
      n_cells = sum(!is.na(GFP_fluorescence) & !is.na(GFP_raw_counts)),
      rho_raw = if (n_cells >= min_n) cor(GFP_fluorescence, GFP_raw_counts, method = "spearman", use = "complete.obs") else NA_real_,
      rho_norm = if (n_cells >= min_n) cor(GFP_fluorescence, GFP_normalized, method = "spearman", use = "complete.obs") else NA_real_,
      reliable = n_cells >= min_n,
      .groups = "drop"
    )
}

overall_rho_raw <- cor(correlation_data$GFP_fluorescence, correlation_data$GFP_raw_counts,
                        method = "spearman", use = "complete.obs")
overall_rho_norm <- cor(correlation_data$GFP_fluorescence, correlation_data$GFP_normalized,
                         method = "spearman", use = "complete.obs")

dox_correlations <- compute_group_correlations(correlation_data, "Dox_Concentration")

message("\nCorrelation Results (Spearman):")
message(sprintf("Overall: raw counts rho = %.3f | normalized (scran) rho = %.3f", overall_rho_raw, overall_rho_norm))
message("\nBy Dox concentration (groups with n < ", MIN_CELLS_FOR_COR, " cells flagged unreliable):")
print(dox_correlations)

if (any(!dox_correlations$reliable)) {
  warning("Some Dox_Concentration groups have too few cells for a reliable correlation estimate -- see 'reliable' column above.")
}

## ---- Scatter plots: raw counts vs. fluorescence, and normalized vs. fluorescence --

# FIX: drop non-positive values before log10 scaling, and log how many were
# dropped, instead of letting ggplot2 silently drop them via the log scale.
make_log_safe <- function(df, cols) {
  n0 <- nrow(df)
  df <- df %>% filter(if_all(all_of(cols), ~ !is.na(.x) & .x > 0))
  n_dropped <- n0 - nrow(df)
  if (n_dropped > 0) {
    message(sprintf("Dropped %d / %d cells with non-positive or missing values before log10 scaling.", n_dropped, n0))
  }
  df
}

plot_gfp_scatter <- function(df, x_col, x_label, title, filename, color_by_dox = TRUE) {
  df_plot <- make_log_safe(df, c(x_col, "GFP_fluorescence"))
  rho <- cor(df_plot[[x_col]], df_plot$GFP_fluorescence, method = "spearman", use = "complete.obs")

  p <- ggplot(df_plot, aes(x = .data[[x_col]], y = GFP_fluorescence))
  if (color_by_dox) {
    p <- p + geom_point(aes(fill = Dox_Concentration, color = Dox_Concentration),
                         alpha = 0.8, size = 2, shape = 21, stroke = 0.2) +
      scale_fill_brewer(palette = "Dark2") + scale_color_brewer(palette = "Dark2")
  } else {
    p <- p + geom_point(alpha = 0.5, size = 2)
  }
  p <- p +
    geom_smooth(method = "lm", se = FALSE, color = "black") +
    scale_x_log10() + scale_y_log10() +
    labs(x = paste0(x_label, " (log10)"), y = "GFP fluorescence (log10)",
         title = paste0(title, "\nOverall Spearman's rho = ", round(rho, 3))) +
    theme(legend.position = "right")
  ggsave(file.path(PLOT_DIR, filename), p, width = 12, height = 8, dpi = 300)
  p
}

p_raw_scatter <- plot_gfp_scatter(correlation_data, "GFP_raw_counts", "GFP raw count",
                                   paste(egfp_gene_id, "raw count vs. fluorescence"),
                                   "gfp_rawcount_vs_fluorescence.png")

p_norm_scatter <- plot_gfp_scatter(correlation_data, "GFP_normalized", "GFP normalized (scran)",
                                    paste(egfp_gene_id, "scran-normalized count vs. fluorescence"),
                                    "gfp_normcount_vs_fluorescence.png")

## ---- Dox-induced-only scatter with per-group + combined rho annotated -----
# FIX: annotation styling (position, bold/color for the summary row) is now
# derived from the actual group table by name, not a hardcoded position
# vector -- robust to a group being empty or reordered.

induced_data <- correlation_data %>% filter(Dox_Concentration %in% INDUCED_LEVELS)
induced_data_log <- make_log_safe(induced_data, c("GFP_normalized", "GFP_fluorescence"))

induced_group_rho <- induced_data_log %>%
  group_by(Dox_Concentration) %>%
  summarize(n_cells = n(),
            rho = if (n() >= MIN_CELLS_FOR_COR) cor(GFP_normalized, GFP_fluorescence, method = "spearman", use = "complete.obs") else NA_real_,
            .groups = "drop") %>%
  mutate(label_group = as.character(Dox_Concentration), is_summary = FALSE)

combined_rho <- cor(induced_data_log$GFP_normalized, induced_data_log$GFP_fluorescence,
                     method = "spearman", use = "complete.obs")

annotation_table <- bind_rows(
  induced_group_rho,
  data.frame(label_group = "All Dox-treated", rho = combined_rho, n_cells = nrow(induced_data_log), is_summary = TRUE)
) %>%
  filter(!is.na(rho)) %>%
  arrange(is_summary) %>%                      # summary row last -> drawn at bottom
  mutate(
    label = sprintf("%s: n=%d, rho = %.3f", label_group, n_cells, rho),
    vjust_pos = 1.5 + (row_number() - 1),        # dynamic instead of hardcoded seq()
    face = if_else(is_summary, "bold", "plain")
  )
# Split into group-level vs. summary rows so each can be drawn with its own
# FIXED (non-aes-mapped) text color, rather than mapping color via aes() and
# adding a second "colour" scale (scale_color_identity()) that would
# conflict with -- and silently override -- the point layer's discrete
# Dox_Concentration color scale for the whole plot.
annotation_groups <- annotation_table %>% filter(!is_summary)
annotation_summary <- annotation_table %>% filter(is_summary)

dox_scatter <- ggplot(induced_data_log, aes(x = GFP_normalized, y = GFP_fluorescence,
                                             fill = Dox_Concentration, color = Dox_Concentration)) +
  geom_point(alpha = 0.8, size = 2, shape = 21, stroke = 0.2) +
  geom_smooth(method = "lm", se = FALSE, color = "black") +
  scale_fill_brewer(palette = "Dark2") + scale_color_brewer(palette = "Dark2") +
  scale_x_log10() + scale_y_log10() +
  labs(x = "GFP normalized (scran, log10)", y = "GFP fluorescence (log10)",
       title = "Normalized GFP vs. fluorescence in Dox-treated cells") +
  theme(legend.position = "right") +
  geom_text(data = annotation_groups,
            aes(label = label, x = Inf, y = Inf, vjust = vjust_pos, fontface = face),
            hjust = 1.1, size = 4, color = "black", inherit.aes = FALSE) +
  geom_text(data = annotation_summary,
            aes(label = label, x = Inf, y = Inf, vjust = vjust_pos, fontface = face),
            hjust = 1.1, size = 4, color = "red", inherit.aes = FALSE)
ggsave(file.path(PLOT_DIR, "gfp_normcount_vs_fluorescence_dox_treated.png"), dox_scatter, width = 12, height = 8, dpi = 300)

## ---- Faceted scatter by Dox concentration (normalized GFP) ----------------
faceted_scatter <- ggplot(induced_data_log, aes(x = GFP_normalized, y = GFP_fluorescence)) +
  geom_point(alpha = 0.5, size = 2) +
  geom_smooth(method = "lm", se = FALSE) +
  facet_wrap(~ Dox_Concentration, scales = "free") +
  scale_x_log10() + scale_y_log10() +
  labs(x = "GFP normalized (scran, log10)", y = "GFP fluorescence (log10)",
       title = "Normalized GFP vs. fluorescence by Dox concentration") +
  geom_text(data = induced_group_rho %>% filter(!is.na(rho)),
            aes(label = sprintf("n=%d\nrho = %.3f", n_cells, rho)),
            x = Inf, y = Inf, hjust = 1.1, vjust = 1.5, size = 4, inherit.aes = FALSE)
ggsave(file.path(PLOT_DIR, "gfp_normcount_vs_fluorescence_by_dox.png"), faceted_scatter, width = 16, height = 6, dpi = 300)

## ============================================================================
## 4. Save results
## ============================================================================
write.csv(correlation_data, file.path(PROCESSED_DATA_DIR, "gfp_correlation_data.csv"), row.names = FALSE)
saveRDS(cell_metadata, file.path(PROCESSED_DATA_DIR, "cell_metadata.rds"))  # now includes GFP_fluorescence

message("GFP fluorescence QC complete. Correlation table and plots saved to: ", PLOT_DIR)
