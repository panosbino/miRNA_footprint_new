# ==============================================================================
# PSCSR benchmark, step 8: inferred activity vs measured miRNA expression.
#
# Ground truth: unit_expression.rds (normalised miRNA counts summed over the
# unit's members). Each cell line is evaluated ONLY on its own top-N units
# (selection.tsv), inside that line. Per method x set type x (line, unit):
#   rho          Spearman(activity, measured expression) in the line's cells.
#   specificity  rho minus the median Spearman(activity, measured expression of
#                another unit on the SAME line's list). Units whose measured
#                expression co-varies with this unit's in that line
#                (|rho| >= SPEC_MAX_TRUTH_RHO, e.g. miR-17~92 members) are left
#                out of the reference; n_reference says how many remain.
#                Near 0 = the score tracks miRNA expression in general, not
#                its own miRNA.
# Across-line correlations are not computed: between lines, any gene set that
# differs between the lines tracks any line-specific miRNA.
# Tools without target lists: bayesReact (step 5), miReact motif (step 6),
# miDGD (sum of predicted proportions of the unit's members in its model).
#
# Outputs (analysis/PSCSR/benchmark/): benchmark_per_line_unit_pscsr.tsv,
#   benchmark_summary_pscsr.tsv, benchmark_pscsr_heatmaps.pdf, benchmark_pscsr_heatmap_per_line.pdf
# Env: MIDGD_PREFIX (default <OUT_DIR>/midgd/4a)
# ==============================================================================
suppressPackageStartupMessages(library(tidyverse))
source(file.path(Sys.getenv("MIRNA_BASE_DIR", "/cfs/klemming/projects/supr/naiss2024-6-235/miRNA_footprint_new"),
                 "scripts/PSCSR/benchmark/config_pscsr.R"))
MIDGD_PREFIX <- Sys.getenv("MIDGD_PREFIX", file.path(OUT_DIR, "midgd/4a"))
NO_LIST <- "no list"

cells <- read_tsv(file.path(DATA_DIR, "cells.tsv"), show_col_types = FALSE)
units <- read_tsv(file.path(DATA_DIR, "units.tsv"), show_col_types = FALSE)
truth <- readRDS(file.path(DATA_DIR, "unit_expression.rds"))[, cells$cell, drop = FALSE]

read_if <- function(path, reader) { if (file.exists(path)) return(reader(path)); cat("SKIPPED (not found):", path, "\n"); NULL }

# ---- List-based scores ---------------------------------------------------------------------------
lst <- bind_rows(read_if(file.path(OUT_DIR, "scores_ours.rds"), readRDS),
                 read_if(file.path(OUT_DIR, "scores_mireact.rds"), readRDS),
                 read_if(file.path(OUT_DIR, "scores_mitea.csv"), \(f) read_csv(f, show_col_types = FALSE)))
lst <- bind_cols(lst, split_key(lst$target_set)) |> dplyr::select(method, unit, set_type, cell, score)

# ---- No-list tools --------------------------------------------------------------------------------
nl <- list()
for (u in units$unit) {
  b <- read_if(file.path(OUT_DIR, "bayesreact", paste0(u, ".rds")), readRDS)
  if (!is.null(b)) {
    nl[[length(nl) + 1]] <- tibble(method = "bayesReact", unit = u, cell = names(b$activity), score = b$activity)
    if (max(b$diagnostics$Rhat, na.rm = TRUE) > 1.01) cat(sprintf("*** bayesReact %s: max Rhat %.3f\n", u, max(b$diagnostics$Rhat, na.rm = TRUE)))
  }
}
trk <- read_if(file.path(OUT_DIR, "mireact_motif_job_tracking.rds"), readRDS)
if (!is.null(trk)) {
  ma <- read_if(file.path(trk$wd, trk$out_file), readRDS)
  if (!is.null(ma)) for (i in seq_len(nrow(units))) {
    if (!units$motif[i] %in% rownames(ma)) { cat("*** motif not in miReact output:", units$motif[i], "\n"); next }
    nl[[length(nl) + 1]] <- tibble(method = "miReact (motif)", unit = units$unit[i], cell = colnames(ma), score = ma[units$motif[i], ])
  }
}
mg <- read_if(paste0(MIDGD_PREFIX, "_mirna_prop.tsv"), \(f) read_tsv(f, show_col_types = FALSE) |> rename(cell = 1))
if (!is.null(mg)) for (i in seq_len(nrow(units))) {
  mem <- intersect(strsplit(units$members[i], ";")[[1]], names(mg))
  if (!length(mem)) { cat("*** miDGD model has no miRNA of unit", units$unit[i], "\n"); next }
  nl[[length(nl) + 1]] <- tibble(method = "miDGD (4a)", unit = units$unit[i], cell = mg$cell,
                                 score = rowSums(as.matrix(mg[, mem, drop = FALSE])))
}
scores <- bind_rows(lst, bind_rows(nl) |> mutate(set_type = NO_LIST)) |>
  inner_join(cells, by = "cell")

# ---- Evaluation pairs: each cell line x the units on ITS top-N list --------------------------------------
selection <- read_tsv(file.path(DATA_DIR, "selection.tsv"), show_col_types = FALSE) |> dplyr::select(cell_type, rank, unit)
sp <- function(x, y) suppressWarnings(cor(x, y, method = "spearman"))
line_cells <- split(cells$cell, cells$cell_type)
truth_rho <- map(line_cells, function(cc) cor(t(truth[, cc, drop = FALSE]), method = "spearman"))   # measured unit-unit, per line

score_mat <- function(d) { s <- d$score; names(s) <- d$cell; s }
res <- scores |>
  filter(is.finite(score)) |>
  group_by(method, set_type, unit) |>
  group_modify(function(d, k) {
    s <- score_mat(d)
    map_dfr(selection$cell_type[selection$unit == k$unit], function(ct) {
      cc <- intersect(line_cells[[ct]], names(s))
      own <- sp(s[cc], truth[k$unit, cc])
      # reference: the other units on this line's list whose measured expression does not co-vary with this unit's
      lst <- selection$unit[selection$cell_type == ct]
      ref <- setdiff(lst[abs(truth_rho[[ct]][k$unit, lst]) < SPEC_MAX_TRUTH_RHO], k$unit)
      others <- vapply(ref, \(j) sp(s[cc], truth[j, cc]), numeric(1))
      tibble(cell_type = ct, n = length(cc), rho = own, n_reference = length(ref),
             rho_others_median = if (length(ref)) median(others) else NA_real_,
             specificity = if (length(ref)) own - median(others) else NA_real_,
             frac_others_below = if (length(ref)) mean(others < own) else NA_real_)
    })
  }) |> ungroup() |>
  left_join(selection, by = c("cell_type", "unit"))
short <- res |> filter(n < lengths(line_cells)[cell_type])
if (nrow(short)) { cat("*** Rows scored on fewer cells than the line has (non-finite scores dropped):\n"); print(short, n = 20) }

sizes <- read_tsv(file.path(OUT_DIR, "target_sets_summary.tsv"), show_col_types = FALSE) |> dplyr::select(unit, set_type, n_genes)
res <- left_join(res, sizes, by = c("unit", "set_type")) |> arrange(method, set_type, cell_type, rank)
write_tsv(res, file.path(OUT_DIR, "benchmark_per_line_unit_pscsr.tsv"))

ranked <- res |> filter(set_type == NO_LIST | n_genes >= 10)          # sets < 10 genes are not ranked
summary_tbl <- ranked |>
  group_by(method, set_type) |>
  summarise(n_pairs = n(), median_rho = median(rho, na.rm = TRUE), median_specificity = median(specificity, na.rm = TRUE),
            .groups = "drop") |>
  left_join(ranked |> group_by(method, set_type, cell_type) |> summarise(m = median(rho, na.rm = TRUE), .groups = "drop") |>
              pivot_wider(names_from = cell_type, values_from = m, names_prefix = "median_rho_"),
            by = c("method", "set_type")) |>
  arrange(desc(median_rho))
write_tsv(summary_tbl, file.path(OUT_DIR, "benchmark_summary_pscsr.tsv"))
cat(sprintf("\nEvaluation: %d cell-line x unit pairs (top %d per line). Medians over pairs with >= 10 genes in the set:\n",
            nrow(selection), TOP_N_PER_TYPE))
print(summary_tbl, n = 40, width = Inf)

# ---- Heatmaps -------------------------------------------------------------------------------------------------
method_levels <- c(paste0("Ours-", c("A", "B", "C", "C2", "D")), "miReact", "miTEA-HiRes",
                   "miReact (motif)", "bayesReact", "miDGD (4a)")
set_levels <- c(unique(sizes$set_type), NO_LIST)
hm <- function(d, col, title, lims, facet = FALSE) {
  p <- d |>
    mutate(method = factor(method, rev(method_levels)), set_type = factor(set_type, set_levels)) |>
    ggplot(aes(set_type, method, fill = .data[[col]])) +
    geom_tile(colour = "white") +
    geom_text(aes(label = sprintf("%.2f", .data[[col]])), size = 2.4) +
    scale_fill_gradient2(low = "#b2182b", mid = "white", high = "#2166ac", midpoint = 0, limits = lims, oob = scales::squish) +
    theme_bw(base_size = 11) + theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(x = NULL, y = NULL, fill = NULL, title = title,
         subtitle = "Each cell line evaluated on its own top-10 units; sets with < 10 genes excluded")
  if (facet) p + facet_wrap(~ cell_type, ncol = 1) else p
}
per_line <- ranked |> group_by(method, set_type, cell_type) |> summarise(median_rho = median(rho, na.rm = TRUE), .groups = "drop")
pdf(file.path(OUT_DIR, "benchmark_pscsr_heatmaps.pdf"), width = 15, height = 6)
print(hm(summary_tbl, "median_rho", "Activity vs measured miRNA within cell lines: median Spearman over line x unit pairs", c(-0.5, 0.5)))
print(hm(summary_tbl, "median_specificity", "Specificity: own miRNA minus median of the line's other (non-co-varying) miRNAs", c(-0.3, 0.3)))
invisible(dev.off())
ggsave(file.path(OUT_DIR, "benchmark_pscsr_heatmap_per_line.pdf"),
       hm(per_line, "median_rho", "Median Spearman per cell line", c(-0.5, 0.5), facet = TRUE), width = 15, height = 16)
cat(sprintf("Wrote benchmark_per_line_unit_pscsr.tsv, benchmark_summary_pscsr.tsv and heatmaps to %s\n", OUT_DIR))
