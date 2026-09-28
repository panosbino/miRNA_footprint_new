#!/usr/bin/env Rscript

# ==============================================================================
# miR-124 activity vs. GFP correlation, across all combinations of:
#   - target-set cutoff (top N TargetScan-predicted targets by cumulative
#     weighted context++ score): 100, 200, 500, 1000
#   - GFP readout: fluorescence (FACS), scran-normalized counts, raw counts
# -> 4 x 3 = 12 correlation plots (x = log10(GFP + 1), the only x-axis scale
#    used -- an unlogged version was tried and dropped as uninformative for
#    this data), individually + as one combined grid, plus a single CSV with
#    both Spearman and Pearson correlation for every combo.
#
# Ported from Basic_miR124_analysis_combined.Rmd, which only computed this
# for a single hardcoded cutoff (top 200) and had to be copy-pasted 3x (once
# per GFP type) to get even that. This version generalizes both.
#
# FIXES APPLIED vs. the Rmd:
#   1. The Rmd's p1/p2 plot titles reference `cor[3]$p.value` -- `cor` is
#      never assigned (only cor_fluo/cor_raw/cor_norm exist); `cor` is the
#      base R function of that name, and `fn[3]` isn't valid, so this would
#      error the moment those plots were built. Each plot below correctly
#      uses its own cor.test() result.
#   2. The Rmd's norm-GFP quantile filter references `GFP_counts_norm$GFP_norm`,
#      but that data frame's column was actually named `GFP_normalized` two
#      lines above -- a typo that would error on a non-existent column. Fixed
#      by building all three GFP data frames directly (see fix 5) with
#      column names used consistently throughout.
#   3. calculate_activity() now reports how many of the requested target IDs
#      actually matched the counts matrix row names (and warns if the match
#      rate looks low), instead of silently dropping unmatched targets via
#      na.omit() with no visibility into how many were lost.
#   4. Cells with zero total expression across all matched target genes
#      (which would make activity = -log2(0) = Inf) are now detected,
#      warned about, and set to NA (excluded from correlation) rather than
#      silently producing an infinite activity score.
#   5. GFP data frames are now built directly with data.frame(), instead of
#      pull() |> as.data.frame() |> rename(<literal deparsed name>) --
#      simpler and not dependent on dplyr's auto-generated column-name text
#      matching exactly across dplyr versions.
#   6. The 99th-percentile GFP outlier cutoff is now computed on the
#      activity-matched subset (post-merge), not the full metadata column
#      pre-merge -- a minor numerical difference from the Rmd, flagged here
#      since it can shift the cutoff slightly.
#   7. Paths pulled into one config block (the Rmd hardcoded three different
#      absolute paths across three different machines' home directories).
# ==============================================================================

## ---- Config -----------------------------------------------------------------
PROJECT_ROOT  <- "~/Desktop/Projects/miRNA_footprint_new"
COUNTS_PATH   <- file.path(PROJECT_ROOT, "datasets", "HEK_SS3", "processed", "scran_normalized_linear.rds")
METADATA_PATH <- file.path(PROJECT_ROOT, "datasets", "HEK_SS3", "processed", "gfp_correlation_data.csv")

# NOT CONFIRMED: no miR-124 target file was visible in resources/ at the time
# this path was set (only miR-199 and mESC target sets were present). Your
# other resource files use "Targets_combined_..." (single underscore) --
# matched here -- rather than the double-underscore "Targets__combined_..."
# used in the source Rmd this script was ported from, which looked like a
# one-off typo rather than your project's actual naming convention. Update
# this to the real filename once the miR-124 target set is in place.
TARGETS_PATH  <- file.path(PROJECT_ROOT, "resources/human", "Targets__combined_124-3p_124-3p.2_506-3p.rds")

# New subfolder, not yet in the project tree -- created on first run.
OUTPUT_DIR    <- file.path(PROJECT_ROOT, "analysis", "HEK_SS3", "mir124_activity_gfp")
PLOT_DIR      <- file.path(OUTPUT_DIR, "plots")
dir.create(PLOT_DIR, showWarnings = FALSE, recursive = TRUE)

CUTOFFS <- c(100, 200, 500, 1000)
GFP_OUTLIER_QUANTILE <- 0.99
MIN_CELLS_FOR_COR <- 10
MIN_TARGET_MATCH_FRAC <- 0.5

## ---- Libraries ----------------------------------------------------------------
library(tidyverse)
library(patchwork)

## ---- Load data (once) ------------------------------------------------------
counts   <- readRDS(COUNTS_PATH)              # genes x cells matrix (scran-normalized, linear scale)
targets  <- readRDS(TARGETS_PATH)             # TargetScan predicted targets for the miR-124 family
metadata <- read.delim(METADATA_PATH, sep = ",")

if (!is.matrix(counts)) counts <- as.matrix(counts)

# Build clean, directly-named GFP lookup tables (fix 5).
gfp_readouts <- list(
  fluo = list(
    df    = data.frame(GFP = metadata$GFP_fluorescence, row.names = metadata$cell_id),
    label = "GFP fluorescence",
    color = "hotpink"
  ),
  norm = list(
    df    = data.frame(GFP = metadata$GFP_normalized, row.names = metadata$cell_id),
    label = "GFP normalized (scran)",
    color = "darkgreen"
  ),
  raw = list(
    df    = data.frame(GFP = metadata$GFP_raw_counts, row.names = metadata$cell_id),
    label = "GFP raw counts",
    color = "cornflowerblue"
  )
)

## ============================================================================
## Activity score: -log2(target-set expression / its mean across cells) + 1
## ============================================================================
# Strong miRNA activity -> lower target gene expression -> higher activity score.
calculate_activity <- function(counts, targets, cutoff_label = "") {
  matched_targets <- intersect(targets$ensembl_gene_id, rownames(counts))
  match_frac <- length(matched_targets) / nrow(targets)
  message(sprintf("  [%s] %d / %d targets matched counts matrix rownames (%.1f%%)",
                   cutoff_label, length(matched_targets), nrow(targets), 100 * match_frac))
  if (match_frac < MIN_TARGET_MATCH_FRAC) {
    warning(sprintf(
      "[%s] Less than %.0f%% of the target set matched the counts matrix -- check gene ID formatting (e.g. Ensembl version suffixes).",
      cutoff_label, 100 * MIN_TARGET_MATCH_FRAC
    ))
  }
  if (length(matched_targets) == 0) {
    stop(sprintf("[%s] No targets matched the counts matrix at all.", cutoff_label))
  }

  cell_sums <- colSums(counts[matched_targets, , drop = FALSE])

  n_zero <- sum(cell_sums == 0)
  if (n_zero > 0) {
    warning(sprintf(
      "[%s] %d cells have zero total expression across matched targets (activity undefined, -log2(0)) -- setting to NA and excluding from correlation.",
      cutoff_label, n_zero
    ))
    cell_sums[cell_sums == 0] <- NA
  }

  cell_sums_norm <- cell_sums / mean(cell_sums, na.rm = TRUE)
  activity <- -log2(cell_sums_norm) + 1
  data.frame(activity = activity, row.names = names(cell_sums))
}

## ============================================================================
## Correlate activity vs. one GFP readout (data + stats only, no plot)
## ============================================================================
compute_correlation <- function(activity_df, gfp_df, gfp_label, cutoff) {
  merged <- merge(activity_df, gfp_df, by = 0)
  rownames(merged) <- merged$Row.names
  merged <- merged[, c("activity", "GFP")]
  merged <- merged[!is.na(merged$activity) & !is.na(merged$GFP), ]

  # 99th-percentile GFP outlier removal, computed on the matched subset (fix 6).
  gfp_cutoff <- quantile(merged$GFP, GFP_OUTLIER_QUANTILE, na.rm = TRUE)
  merged <- merged[merged$GFP < gfp_cutoff, ]

  n_cells <- nrow(merged)
  if (n_cells < MIN_CELLS_FOR_COR) {
    warning(sprintf(
      "Only %d cells available for %s / top%d targets -- correlation may be unreliable (min recommended: %d).",
      n_cells, gfp_label, cutoff, MIN_CELLS_FOR_COR
    ))
  }

  test <- suppressWarnings(cor.test(merged$activity, merged$GFP, method = "spearman", exact = FALSE))
  spearman_rho <- unname(test$estimate)
  spearman_p   <- test$p.value

  list(
    merged = merged,
    stats = data.frame(gfp_type = gfp_label, cutoff = cutoff,
                        spearman_rho = spearman_rho, spearman_p = spearman_p,
                        n_cells = n_cells)
  )
}

## ============================================================================
## Build the scatter plot from a precomputed correlation result, and compute
## Pearson correlation on the same log10(GFP + 1) scale being plotted.
## ============================================================================
build_plot <- function(merged, stats, gfp_label, color, cutoff) {
  x_values <- log10(merged$GFP + 1)

  pearson_test <- suppressWarnings(cor.test(x_values, merged$activity, method = "pearson"))
  pearson_r <- unname(pearson_test$estimate)
  pearson_p <- pearson_test$p.value

  p <- ggplot(merged, aes(x = log10(GFP + 1), y = activity)) +
    geom_point(color = color, alpha = 0.6) +
    labs(
      title = paste0(gfp_label, " | top ", cutoff, " targets"),
      subtitle = sprintf("Spearman rho = %.3f, p = %s   |   Pearson r = %.3f, p = %s   (n = %d)",
                          stats$spearman_rho, formatC(stats$spearman_p, format = "e", digits = 2),
                          pearson_r, formatC(pearson_p, format = "e", digits = 2),
                          stats$n_cells),
      x = paste0("log10(", gfp_label, " + 1)"),
      y = "miR-124 activity score"
    ) +
    theme_bw(base_size = 13) +
    theme(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5, size = 10))

  list(
    plot = p,
    pearson_stats = data.frame(
      gfp_type = gfp_label, cutoff = cutoff,
      pearson_r = pearson_r, pearson_p = pearson_p, n_cells = stats$n_cells
    )
  )
}

## ============================================================================
## Run all 12 combinations (target-set cutoff x GFP readout)
## ============================================================================
targets_sorted <- targets[order(targets$Cumulative.weighted.context...score), ]

all_plots <- list()
all_spearman_stats <- list()
all_pearson_stats  <- list()

for (cutoff in CUTOFFS) {
  message(sprintf("\nComputing activity for top %d targets...", cutoff))
  if (cutoff > nrow(targets_sorted)) {
    warning(sprintf("Requested cutoff %d exceeds available targets (%d) -- using all available targets instead.",
                     cutoff, nrow(targets_sorted)))
  }
  target_subset <- targets_sorted[seq_len(min(cutoff, nrow(targets_sorted))), ]
  activity <- calculate_activity(counts, target_subset, cutoff_label = paste0("top", cutoff))

  for (gfp_type in names(gfp_readouts)) {
    r <- gfp_readouts[[gfp_type]]
    result <- compute_correlation(activity, r$df, r$label, cutoff)

    key <- paste0("top", cutoff, "_", gfp_type)
    all_spearman_stats[[key]] <- result$stats

    plot_result <- build_plot(result$merged, result$stats, r$label, r$color, cutoff)
    all_plots[[key]] <- plot_result$plot
    all_pearson_stats[[key]] <- plot_result$pearson_stats
    ggsave(
      file.path(PLOT_DIR, paste0("mir124_activity_vs_", gfp_type, "_top", cutoff, ".png")),
      plot_result$plot, width = 7, height = 6, dpi = 300
    )
  }
}

## ---- Combined grid: rows = cutoff, cols = GFP readout ----------------------
combined <- wrap_plots(all_plots, ncol = length(gfp_readouts)) +
  plot_annotation(
    title = sprintf("miR-124 activity vs. GFP readouts across target-set cutoffs (top %.0f%% GFP outliers removed)",
                     100 * GFP_OUTLIER_QUANTILE),
    theme = theme(plot.title = element_text(hjust = 0.5, size = 16))
  )
ggsave(file.path(PLOT_DIR, "mir124_activity_vs_GFP_all_cutoffs_grid.png"),
       combined, width = 18, height = 6 * length(CUTOFFS), dpi = 300, limitsize = FALSE)

combined[1] |> 
  as.data.frame() |>
  ggplot() + 
  geom_point(aes(y = data.activity,x = log10(data.GFP)),fill = "#4cadad" , size = 3, shape = 21, color = "#008a8a" ) +
  labs(
    title = paste0(
      "Log10 GFP vs. Log10 miR-124 counts\n",
      "Spearman \u03c1 = X"
    ),
    x = "log10 GFP",
    y = "miRNA activity score"
  ) +
  theme_bw(base_size = 13) +
  theme(
    plot.title = element_text(hjust = 0.5)
  )

## ---- Summary table of all correlations ---------------------------------------
# One row per (cutoff, gfp_type) for each method -- Spearman and Pearson both
# computed on the log10(GFP + 1) scale actually plotted.
spearman_df <- bind_rows(all_spearman_stats) %>%
  mutate(method = "spearman", estimate = spearman_rho, p_value = spearman_p) %>%
  dplyr::select(gfp_type, cutoff, method, estimate, p_value, n_cells)

pearson_df <- bind_rows(all_pearson_stats) %>%
  mutate(method = "pearson", estimate = pearson_r, p_value = pearson_p) %>%
  dplyr::select(gfp_type, cutoff, method, estimate, p_value, n_cells)

stats_df <- bind_rows(spearman_df, pearson_df) %>%
  arrange(cutoff, gfp_type, method)

print(stats_df)
write.csv(stats_df, file.path(OUTPUT_DIR, "mir124_activity_gfp_correlations.csv"), row.names = FALSE)

message("\nDone. 12 individual plots + 1 combined grid + Spearman+Pearson correlation summary CSV saved to: ", OUTPUT_DIR)

