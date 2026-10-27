# ==============================================================================
# PSCSR diagnostic: why does activity computed from LOG-normalised expression
# (sum of log1p values, as in the old Analysis.Rmd) correlate more with measured
# miRNA expression than activity from LINEAR normalised expression?
#
# Hypothesis: after normalisation every cell has the same total. A linear sum
# over a large gene set is then close to a fixed fraction of that total, but a
# sum of logs is not constrained: for a fixed total, sum(log1p(x)) grows the more
# evenly the total is spread over genes (log is concave). Log-sum activity would
# then track transcriptome complexity (genes detected) for ANY gene set, and
# correlate with miRNA expression whenever miRNA expression also tracks it.
#
# Parts (all within cell lines, on the benchmark's (line, unit) pairs):
#   1. Cell covariates vs activity of RANDOM gene sets, log-sum vs linear-sum.
#   2. Cell covariates vs measured miRNA expression.
#   3. Activity vs miRNA, targets and random sets, before and after removing
#      each covariate from both sides (partial Spearman).
#   4. Direct test of the concavity explanation: log-sum activity of random sets
#      vs an evenness index of the cell's transcriptome.
#
# Needs steps 0 and 1 of the benchmark (data/ and target_sets.rds).
# Env: MIRNA_BASE_DIR; SET_TYPE (default TargetScan_all); N_RANDOM (default 20)
# Output: analysis/PSCSR/benchmark/diagnostics/log_vs_linear_*.tsv
# ==============================================================================
suppressPackageStartupMessages({ library(tidyverse); library(data.table) })
source(file.path(Sys.getenv("MIRNA_BASE_DIR", "/cfs/klemming/projects/supr/naiss2024-6-235/miRNA_footprint_new"),
                 "scripts/PSCSR/benchmark/config_pscsr.R"))
SET_TYPE <- Sys.getenv("SET_TYPE", "TargetScan_all")
N_RANDOM <- as.integer(Sys.getenv("N_RANDOM", "20"))
DIAG_DIR <- file.path(OUT_DIR, "diagnostics"); dir.create(DIAG_DIR, showWarnings = FALSE)
set.seed(20261021)

cells <- read_tsv(file.path(DATA_DIR, "cells.tsv"), show_col_types = FALSE)
sel   <- read_tsv(file.path(DATA_DIR, "selection.tsv"), show_col_types = FALSE)
truth <- readRDS(file.path(DATA_DIR, "unit_expression.rds"))[, cells$cell]
ts    <- readRDS(file.path(OUT_DIR, "target_sets.rds"))
U     <- ts$universe

mr <- load_mrna_counts()[, cells$cell]; mr <- mr[sub("\\.\\d+$", "", rownames(mr)) %in% U, ]
rownames(mr) <- sub("\\.\\d+$", "", rownames(mr))
mi <- load_mirna_counts()[, cells$cell]
lib_all <- colSums(load_mrna_counts()[, cells$cell])                     # all genes, as Seurat would use

# Same genes (universe U), two transforms of the SAME normalised values
LIN <- sweep(mr, 2, lib_all, "/") * 1e4                                 # linear, per-cell scaled to 10,000
LOG <- log1p(LIN)                                                       # = Seurat LogNormalize (GEO tables)
act <- function(M, g) { s <- colSums(M[g, , drop = FALSE]); -log2(s / mean(s)) + 1 }   # Utils.R formula

# Cell covariates
p <- LIN / colSums(LIN)
cov <- tibble(cell = cells$cell, cell_type = cells$cell_type,
              mrna_umis = lib_all, genes_detected = colSums(mr > 0),
              mirna_umis = colSums(mi), mirnas_detected = colSums(mi > 0),
              evenness = -colSums(ifelse(p > 0, p * log(p), 0)) / log(colSums(mr > 0)))   # Pielou evenness of the transcriptome
COVS <- c("mrna_umis", "genes_detected", "mirna_umis", "mirnas_detected", "evenness")
sp <- function(x, y) cor(x, y, method = "spearman")
pairs <- sel |> mutate(key = set_key(unit, SET_TYPE), targets = map(key, \(k) ts$sets[[k]])) |>
  filter(lengths(targets) >= 10)
cat(sprintf("Set type %s: %d of %d (line, unit) pairs have >= 10 targets; median %d targets\n",
            SET_TYPE, nrow(pairs), nrow(sel), as.integer(median(lengths(pairs$targets)))))

# ---- Part 1: random gene sets vs cell covariates ------------------------------------------------
rand_sets <- replicate(N_RANDOM, sample(U, as.integer(median(lengths(pairs$targets)))), simplify = FALSE)
part1 <- map_dfr(c("log_sum" = "LOG", "linear_sum" = "LIN"), .id = "activity", function(nm) {
  M <- get(nm)
  map_dfr(split(seq_len(nrow(cov)), cov$cell_type), .id = "cell_type", function(i)
    map_dfr(rand_sets, function(g) { a <- act(M, g)[i]; as_tibble(as.list(sapply(COVS, \(cv) sp(a, cov[[cv]][i])))) }) |>
      summarise(across(everything(), median)))
})
cat("\n== Part 1: random gene sets -- Spearman(activity, cell covariate) within lines (median over sets) ==\n")
print(part1 |> mutate(across(where(is.numeric), \(x) round(x, 2))), n = Inf, width = Inf)

# ---- Part 2: measured miRNA vs cell covariates -------------------------------------------------------
part2 <- pairs |> rowwise() |>
  mutate(res = list({ i <- cells$cell_type == cell_type; as_tibble(as.list(sapply(COVS, \(cv) sp(truth[unit, i], cov[[cv]][i])))) })) |>
  ungroup() |> dplyr::select(cell_type, unit, res) |> unnest(res)
cat("\n== Part 2: Spearman(measured miRNA unit, cell covariate) within lines, median over pairs per line ==\n")
print(part2 |> group_by(cell_type) |> summarise(across(all_of(COVS), median)) |> mutate(across(where(is.numeric), \(x) round(x, 2))), width = Inf)

# ---- Part 3: activity vs miRNA, raw and partial -------------------------------------------------------
rank_resid <- function(y, z) resid(lm(rank(y) ~ rank(z)))               # partial Spearman via rank residuals
part3 <- pmap_dfr(pairs, function(cell_type, unit, targets, ...) {
  i <- cells$cell_type == cell_type; y <- truth[unit, i]
  rs <- replicate(5, sample(U, length(targets)), simplify = FALSE)
  map_dfr(c(log_sum = "LOG", linear_sum = "LIN"), .id = "activity", function(nm) {
    M <- get(nm); at <- act(M, targets)[i]; ar <- lapply(rs, \(g) act(M, g)[i])
    map_dfr(c("none", COVS), function(cv) {
      f <- if (cv == "none") \(v) v else \(v) rank_resid(v, cov[[cv]][i])
      tibble(controlled_for = cv, targets = sp(f(at), f(y)),
             random = median(sapply(ar, \(a) sp(f(a), f(y)))))
    })
  }) |> mutate(cell_type = cell_type, unit = unit)
})
part3_sum <- part3 |> group_by(activity, controlled_for) |>
  summarise(targets = median(targets), random = median(random), gap = median(targets - random), .groups = "drop") |>
  arrange(activity, factor(controlled_for, c("none", COVS)))
cat("\n== Part 3: Spearman(activity, measured miRNA), median over pairs; covariate removed from both sides ==\n")
print(part3_sum |> mutate(across(where(is.numeric), \(x) round(x, 2))), n = Inf)

# ---- Part 4: concavity check -------------------------------------------------------------------------------
cat("\n== Part 4: log-sum activity of random sets vs transcriptome evenness and genes detected, within lines ==\n")
cat("(negative = more even / more complex transcriptome -> larger sum of logs -> lower activity)\n")
print(part1 |> filter(activity == "log_sum") |> dplyr::select(cell_type, evenness, genes_detected) |>
        mutate(across(where(is.numeric), \(x) round(x, 2))))

write_tsv(part1, file.path(DIAG_DIR, "log_vs_linear_part1_random_vs_covariates.tsv"))
write_tsv(part2, file.path(DIAG_DIR, "log_vs_linear_part2_mirna_vs_covariates.tsv"))
write_tsv(part3, file.path(DIAG_DIR, "log_vs_linear_part3_partial.tsv"))
cat(sprintf("\nWrote tables to %s\n", DIAG_DIR))
