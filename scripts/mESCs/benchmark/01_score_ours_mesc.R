# ==============================================================================
# mESC benchmark, step 1: our method, 5 scoring variants x every target set.
#
# Variants, exactly as in "Plan for benchmarking" (higher = more miRNA activity):
#   A  raw sum                       calculate_activity() from scripts/Utils.R
#   B  equal weight per gene         each gene divided by its own mean first
#   C  mean gene-centred log ratio
#   C2 C minus the non-target mean   background = genes in U that are in no
#                                    mESC-family target set of any database
#   D  mean within-cell rank         of the gene-centred log ratios
#
# A here normalises by the mean over ALL cells (Utils.R), whereas
# mESC_our_method_separation.R normalises by the KO mean. The two differ by a
# constant shift, so the separation metric (which only uses the order of cells)
# is identical.
#
# Input : analysis/mESCs/benchmark/target_sets.rds  (step 0)
# Output: analysis/mESCs/benchmark/scores_ours.rds   long table: method, target_set, cell, score
#
# Environment: MIRNA_BASE_DIR, PSEUDOCOUNT (default 1, in CPM units)
# ==============================================================================

suppressPackageStartupMessages(library(tidyverse))

BASE_DIR    <- Sys.getenv("MIRNA_BASE_DIR", "/cfs/klemming/projects/supr/naiss2024-6-235/miRNA_footprint_new")
OUT_DIR     <- file.path(BASE_DIR, "analysis/mESCs/benchmark")
PSEUDOCOUNT <- as.numeric(Sys.getenv("PSEUDOCOUNT", "1"))
setwd(BASE_DIR)
source("scripts/Utils.R")   # calculate_activity()

ts_obj <- readRDS(file.path(OUT_DIR, "target_sets.rds"))
SETS <- ts_obj$sets

pheno <- read.csv("datasets/mESCs/phenotype.csv", fileEncoding = "UTF-8-BOM")
stopifnot("cells" %in% names(pheno))
counts <- as.matrix(readRDS("datasets/mESCs/final_counts.rds"))
X <- counts[ts_obj$universe, colnames(counts) %in% pheno$cells]
cat(sprintf("X: %d genes (universe U) x %d cells\n", nrow(X), ncol(X)))

# ---- Precomputed transforms (same definitions as tarbase_diag_common.R) ---------
gene_mean   <- rowMeans(X)
stopifnot(all(gene_mean > 0))
LOGRATIO    <- log2((X + PSEUDOCOUNT) / (gene_mean + PSEUDOCOUNT))
RANKS       <- apply(LOGRATIO, 2, rank, ties.method = "average")
NON_TARGETS <- setdiff(rownames(X), ts_obj$focal_union)
cat(sprintf("Non-target background for C2: %d genes\n", length(NON_TARGETS)))
non_target_mean <- colMeans(LOGRATIO[NON_TARGETS, , drop = FALSE])

score_variant <- function(genes, variant) {
  genes <- intersect(genes, rownames(X))
  s <- switch(variant,
    A  = calculate_activity(counts = X, targets = genes),
    B  = { v <- colSums(X[genes, , drop = FALSE] / gene_mean[genes]); -log2(v / mean(v)) + 1 },
    C  = -colMeans(LOGRATIO[genes, , drop = FALSE]),
    C2 = -(colMeans(LOGRATIO[genes, , drop = FALSE]) - non_target_mean),
    D  = -colMeans(RANKS[genes, , drop = FALSE]),
    stop("Unknown variant ", variant))
  setNames(as.numeric(s), colnames(X))
}
VARIANTS <- c("A", "B", "C", "C2", "D")

scores <- map_dfr(VARIANTS, function(v) imap_dfr(SETS, function(genes, set_name) {
  s <- score_variant(genes, v)
  tibble(method = paste0("Ours-", v), target_set = set_name, cell = names(s), score = s)
}))

bad <- scores |> group_by(method, target_set) |> summarise(n_nonfinite = sum(!is.finite(score)), .groups = "drop") |> filter(n_nonfinite > 0)
if (nrow(bad) > 0) {
  cat("Non-finite scores (A/B give +Inf when a cell has zero expression of every target; kept, ranked as highest activity):\n")
  print(bad, n = Inf)
}

saveRDS(scores, file.path(OUT_DIR, "scores_ours.rds"))
cat(sprintf("Wrote scores_ours.rds: %d variants x %d target sets x %d cells\n",
            length(VARIANTS), length(SETS), ncol(X)))
