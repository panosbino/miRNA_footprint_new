# ==============================================================================
# HEK293 benchmark, step 5: accuracy of every method x target set against GFP
# fluorescence (FACS index sort, FITC-A).
#
# Ground truth and cells: load_gfp_truth() from scripts/Utils.R, reading
# datasets/HEK_SS3/processed/gfp_correlation_data.csv (fluorescence already
# joined to cells). Cells = its QC-passing set (GFP mRNA below the 99th
# percentile), because the saved bayesReact and miReact-motif tables only cover
# those cells; every method is scored on the same cells.
# Cell sets: "all" QC-passing cells and "induced" (0.01, 0.1, 1 ug/ml Dox).
# Metric: Spearman rho (higher activity should mean higher GFP), approximate
# 95% CI (Bonett & Wright 2000, as in fluorescence_accuracy_comparison.R).
# NA scores (miTEA-HiRes p = 0) are dropped per row and counted in n_dropped.
#
# Tools without target lists, run as-is: bayesReact and miReact motif mode
# (saved results in analysis/HEK_SS3/comparisons/), miDGD (step 4).
#
# Outputs (analysis/HEK_SS3/benchmark/):
#   benchmark_accuracy_hek.tsv, benchmark_accuracy_hek_heatmap.pdf
# ==============================================================================

suppressPackageStartupMessages(library(tidyverse))

BASE_DIR <- Sys.getenv("MIRNA_BASE_DIR", "/cfs/klemming/projects/supr/naiss2024-6-235/miRNA_footprint_new")
PROC_DIR <- file.path(BASE_DIR, "datasets/HEK_SS3/processed")
COMP_DIR <- file.path(BASE_DIR, "analysis/HEK_SS3/comparisons")
OUT_DIR  <- file.path(BASE_DIR, "analysis/HEK_SS3/benchmark")
MIDGD_OUT   <- Sys.getenv("MIDGD_PREFIX", file.path(OUT_DIR, "midgd/4a"))
MIDGD_LABEL <- "miDGD (4a)"
setwd(BASE_DIR)
source("scripts/Utils.R")   # load_gfp_truth()

truth <- load_gfp_truth(PROC_DIR, upper_pct = 0.99)
truth <- truth[truth$qc_pass, ]

# ---- Scores ----------------------------------------------------------------------------
read_if <- function(path, reader) {
  if (file.exists(path)) return(reader(path))
  cat("SKIPPED (not found):", path, "\n"); NULL
}
scores <- bind_rows(
  read_if(file.path(OUT_DIR, "scores_ours.rds"), readRDS),
  read_if(file.path(OUT_DIR, "scores_mireact.rds"), readRDS),
  read_if(file.path(OUT_DIR, "scores_mitea.csv"), \(f) read_csv(f, show_col_types = FALSE))
)

NO_LIST <- "no list (motif / expression)"
saved <- list(bayesReact        = c(file = "res_bayesreact_comparison_scran.rds",       el = "bayesreact"),
              `miReact (motif)` = c(file = "res_mireact_default_comparison_scran.rds",  el = "mireact_default"))
for (nm in names(saved)) {
  r <- read_if(file.path(COMP_DIR, saved[[nm]][["file"]]), readRDS)
  if (is.null(r)) next
  d <- r[[saved[[nm]][["el"]]]]
  col <- setdiff(names(d)[sapply(d, is.numeric)], "eGFP")[1]
  scores <- bind_rows(scores, tibble(method = nm, target_set = NO_LIST, cell = d$Cell_ID, score = d[[col]]))
}
midgd <- read_if(paste0(MIDGD_OUT, "_mirna_prop.tsv"),
                 \(f) read_tsv(f, show_col_types = FALSE) |> rename(cell = 1))
if (!is.null(midgd)) {
  stopifnot("hsa-miR-124-3p" %in% names(midgd))
  scores <- bind_rows(scores, tibble(method = MIDGD_LABEL, target_set = NO_LIST,
                                     cell = midgd$cell, score = midgd$`hsa-miR-124-3p`))
  cat(sprintf("miDGD predicted miR-124-3p proportion: median %.2e, IQR/median %.2f (near-constant predictions make rho noise)\n",
              median(midgd$`hsa-miR-124-3p`), IQR(midgd$`hsa-miR-124-3p`) / median(midgd$`hsa-miR-124-3p`)))
}

# ---- Metric ---------------------------------------------------------------------------------
spearman_ci <- function(r, n) { se <- sqrt((1 + r^2 / 2) / (n - 3)); tanh(atanh(r) + c(-1, 1) * 1.96 * se) }
eval_one <- function(score, fluor) {
  ok <- is.finite(score); n <- sum(ok)
  r <- cor(score[ok], fluor[ok], method = "spearman")
  ci <- spearman_ci(r, n)
  tibble(n = n, n_dropped = sum(!ok), rho = r, ci_low = ci[1], ci_high = ci[2])
}
d <- inner_join(scores, truth, by = c("cell" = "Cell_ID"))
res <- bind_rows(
  d |> mutate(cell_set = "all"),
  d |> filter(induced) |> mutate(cell_set = "induced")
) |>
  group_by(method, target_set, cell_set) |>
  summarise(eval_one(score, gfp_facs), .groups = "drop")

expected_n <- c(all = nrow(truth), induced = sum(truth$induced))
short <- res |> filter(n + n_dropped != expected_n[cell_set])
if (nrow(short) > 0) { cat("*** Rows not covering every evaluation cell:\n"); print(short, n = Inf) }
if (any(res$n_dropped > 0)) { cat("*** Cells dropped (NA score):\n"); print(filter(res, n_dropped > 0), n = Inf) }

ref <- sapply(c(all = "all", induced = "induced"), function(cs) {
  t <- if (cs == "induced") truth[truth$induced, ] else truth
  cor(t$gfp_mrna, t$gfp_facs, method = "spearman")
})
cat(sprintf("\nReference: GFP mRNA vs fluorescence on the same cells: all %.3f, induced %.3f\n", ref[["all"]], ref[["induced"]]))

sizes <- read_tsv(file.path(OUT_DIR, "target_sets_summary.tsv"), show_col_types = FALSE)
res <- left_join(res, dplyr::select(sizes, target_set, n_genes), by = "target_set")
write_tsv(res, file.path(OUT_DIR, "benchmark_accuracy_hek.tsv"))
print(res |> select(-n_dropped) |> pivot_wider(names_from = cell_set, values_from = c(rho, ci_low, ci_high, n)) |>
        arrange(desc(rho_induced)), n = Inf, width = Inf)

# ---- Heatmap ---------------------------------------------------------------------------------
method_levels <- c(paste0("Ours-", c("A", "B", "C", "C2", "D")), "miReact", "miTEA-HiRes",
                   "miReact (motif)", "bayesReact", MIDGD_LABEL)
set_levels <- c(sizes$target_set, NO_LIST)
set_label  <- \(x) ifelse(x %in% sizes$target_set,
                          sprintf("%s (%d)", x, sizes$n_genes[match(x, sizes$target_set)]), x)
p <- res |>
  mutate(method = factor(method, rev(method_levels)), target_set = factor(target_set, set_levels),
         cell_set = factor(cell_set, c("induced", "all"),
                           labels = sprintf("%s cells (n = %d; GFP mRNA vs fluorescence = %.2f)",
                                            c("Dox-induced", "All QC"), expected_n[c("induced", "all")],
                                            ref[c("induced", "all")]))) |>
  ggplot(aes(target_set, method, fill = rho)) +
  geom_tile(colour = "white") +
  geom_text(aes(label = sprintf("%.2f", rho)), size = 2.4) +
  scale_fill_gradient2(low = "#b2182b", mid = "white", high = "#2166ac", midpoint = 0, limits = c(-1, 1),
                       name = "Spearman rho\nvs fluorescence") +
  scale_x_discrete(labels = set_label) +
  facet_wrap(~ cell_set, ncol = 1) +
  theme_bw(base_size = 11) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = NULL, y = NULL, title = "HEK293 miR-124-3p activity vs GFP fluorescence",
       subtitle = "Sets with < 10 genes should not be ranked; differences within the CIs are not a ranking")
ggsave(file.path(OUT_DIR, "benchmark_accuracy_hek_heatmap.pdf"), p, width = 17, height = 9)
cat(sprintf("Wrote benchmark_accuracy_hek.tsv and benchmark_accuracy_hek_heatmap.pdf to %s\n", OUT_DIR))
