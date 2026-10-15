# ==============================================================================
# mESC benchmark, step 4: KO-vs-Control separation for every method x target set.
#
# Metric: sep_error() / sep_error_cv() from scripts/mESCs/mESC_separation_utils.R
# (balanced error at the best single cutoff, Control = higher activity; 5-fold CV).
# Raw scores are used, so the CV error matches the existing pipeline (the CV
# cutoff is a midpoint between training cells, so it is NOT rank-invariant).
# +/-Inf (A/B with zero target expression) is replaced by a value just beyond
# the finite range, which keeps the order of cells.
# NA scores (miTEA-HiRes p = 0) are dropped per (method, target set) and counted
# in n_dropped, which is printed.
#
# Also includes the saved "combined" results of the tools run without target
# lists (bayesReact, miReact motif mode); they are not rerun.
#
# Outputs (analysis/mESCs/benchmark/):
#   benchmark_separation_mesc.tsv, benchmark_separation_mesc_heatmap.pdf
# ==============================================================================

suppressPackageStartupMessages(library(tidyverse))

BASE_DIR <- Sys.getenv("MIRNA_BASE_DIR", "/cfs/klemming/projects/supr/naiss2024-6-235/miRNA_footprint_new")
OUT_DIR  <- file.path(BASE_DIR, "analysis/mESCs/benchmark")
setwd(BASE_DIR)
u <- new.env()   # sourced into its own environment so its hard-coded BASE_DIR does not overwrite ours
sys.source("scripts/mESCs/mESC_separation_utils.R", envir = u)

pheno <- read.csv("datasets/mESCs/phenotype.csv", fileEncoding = "UTF-8-BOM")
stopifnot("cells" %in% names(pheno))
pheno$Exp[pheno$Exp == "WT"] <- "Control"

# ---- Scores from steps 1-3 --------------------------------------------------------
read_if <- function(path, reader) {
  if (file.exists(path)) return(reader(path))
  cat("SKIPPED (not found):", path, "\n"); NULL
}
scores <- bind_rows(
  read_if(file.path(OUT_DIR, "scores_ours.rds"), readRDS),
  read_if(file.path(OUT_DIR, "scores_mireact.rds"), readRDS),
  read_if(file.path(OUT_DIR, "scores_mitea.csv"), \(f) read_csv(f, show_col_types = FALSE))
)

# ---- Tools without target lists: saved "combined" results ---------------------------
AS_IS <- c(bayesReact = "res_bayesreact_separation.rds", `miReact (motif)` = "res_mireact_separation.rds")
for (nm in names(AS_IS)) {
  m <- read_if(file.path(BASE_DIR, "analysis/mESCs", AS_IS[[nm]]), readRDS)$combined$merged
  if (!is.null(m)) scores <- bind_rows(scores, tibble(method = nm, target_set = "motif (no list)",
                                                      cell = m$cells, score = m$activity))
}

# ---- Metric -----------------------------------------------------------------------------
eval_one <- function(score, exp) {
  ok <- !is.na(score)
  r <- score[ok]; is_ko <- exp[ok] == "KO"
  fin <- r[is.finite(r)]; pad <- diff(range(fin)) + 1
  r[r ==  Inf] <- max(fin) + pad
  r[r == -Inf] <- min(fin) - pad
  s <- u$sep_error(r, is_ko)
  tibble(n_cells = sum(ok), n_dropped = sum(!ok), cells_wrong = s$cells_wrong,
         ctrl_called_ko = s$ctrl_called_ko, ko_called_ctrl = s$ko_called_ctrl,
         bal_error = 100 * s$bal_error, bal_error_cv = 100 * u$sep_error_cv(r, is_ko))
}
res <- scores |>
  inner_join(pheno, by = c("cell" = "cells")) |>
  group_by(method, target_set) |>
  summarise(eval_one(score, Exp), .groups = "drop")

if (any(res$n_dropped > 0)) { cat("*** Cells dropped (NA score):\n"); print(filter(res, n_dropped > 0), n = Inf) }
if (any(res$n_cells + res$n_dropped != nrow(pheno)))
  cat("*** WARNING: some (method, target set) rows do not cover all", nrow(pheno), "cells in phenotype.csv\n")

sizes <- read_tsv(file.path(OUT_DIR, "target_sets_summary.tsv"), show_col_types = FALSE)
res <- left_join(res, dplyr::select(sizes, target_set, n_genes), by = "target_set")
write_tsv(res, file.path(OUT_DIR, "benchmark_separation_mesc.tsv"))
print(arrange(res, bal_error), n = Inf, width = Inf)

# ---- Heatmap ------------------------------------------------------------------------------
method_levels <- c(paste0("Ours-", c("A", "B", "C", "C2", "D")), "miReact", "miTEA-HiRes",
                   "miReact (motif)", "bayesReact")
set_levels <- c(sizes$target_set, "motif (no list)")
set_label  <- \(x) ifelse(x %in% sizes$target_set,
                          sprintf("%s (%d)", x, sizes$n_genes[match(x, sizes$target_set)]), x)
p <- res |>
  mutate(method = factor(method, rev(method_levels)), target_set = factor(target_set, set_levels)) |>
  ggplot(aes(target_set, method, fill = bal_error)) +
  geom_tile(colour = "white") +
  geom_text(aes(label = sprintf("%.1f", bal_error)), size = 2.6) +
  scale_fill_viridis_c(direction = -1, limits = c(0, 50), oob = scales::squish, name = "Balanced\nerror (%)") +
  scale_x_discrete(labels = set_label) +
  theme_bw(base_size = 11) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = NULL, y = NULL,
       title = "mESC Drosha-KO vs Control: balanced misclassification at the best cutoff",
       subtitle = "50% = chance; negctrl_TargetScan should be near 50%; sets < 10 genes should not be ranked")
ggsave(file.path(OUT_DIR, "benchmark_separation_mesc_heatmap.pdf"), p, width = 15, height = 5.5)
cat(sprintf("Wrote benchmark_separation_mesc.tsv and benchmark_separation_mesc_heatmap.pdf to %s\n", OUT_DIR))
