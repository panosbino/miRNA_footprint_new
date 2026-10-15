# ==============================================================================
# mESC benchmark, step 2: miReact with each target set.
#
# miReact's TarBase mode (addMotifs() with runparameters$tarbaserun = TRUE) does
# two things: a constant motif p-value of 0.01 for every gene, and a 0/1
# membership vector per miRNA. wcmod.p() is then run per cell on genes ordered
# by (expression - gene median), decreasing. This script runs the same
# wcmod.p() call with the same ordering and the same constant p, swapping in
# each benchmark target set as the membership vector, so any list can be
# scored the way miReact scores TarBase. Sign flipped (* -1) as in wrapper3.R,
# so higher = targets lower = more activity.
#
# Gene universe: U from step 0 (all genes in U have a miReact 3'UTR entry, so
# addSeqs() would keep every one of them; it is not called).
#
# Input : analysis/mESCs/benchmark/target_sets.rds, tools/miReact/code/wcmod.p.R
# Output: analysis/mESCs/benchmark/scores_mireact.rds
# ==============================================================================

suppressPackageStartupMessages(library(tidyverse))

BASE_DIR    <- Sys.getenv("MIRNA_BASE_DIR", "/cfs/klemming/projects/supr/naiss2024-6-235/miRNA_footprint_new")
MIREACT_DIR <- file.path(BASE_DIR, "tools/miReact")
OUT_DIR     <- file.path(BASE_DIR, "analysis/mESCs/benchmark")
CONST_P     <- 0.01     # value addMotifs() puts in pval.mat in TarBase mode
ALPHA       <- 1e-10    # as in mireact_comparison_rebuild.R
setwd(BASE_DIR)
source(file.path(MIREACT_DIR, "code/wcmod.p.R"))

ts_obj <- readRDS(file.path(OUT_DIR, "target_sets.rds"))
SETS <- ts_obj$sets

pheno <- read.csv("datasets/mESCs/phenotype.csv", fileEncoding = "UTF-8-BOM")
stopifnot("cells" %in% names(pheno))
counts <- as.matrix(readRDS("datasets/mESCs/final_counts.rds"))
X <- counts[ts_obj$universe, colnames(counts) %in% pheno$cells]
cat(sprintf("X: %d genes x %d cells\n", nrow(X), ncol(X)))

# ---- miReact's per-cell ordering (as in mireact_comparison_rebuild.R) ----------------
medianExp  <- apply(X, 1, median)
rank_order <- apply(X, 2, function(x) order(x - medianExp, decreasing = TRUE))
p_vec      <- rep(CONST_P, nrow(X))
n_tied <- colSums(sweep(X, 1, medianExp) == 0)
cat(sprintf("Genes exactly at their median (tied; order() keeps row order among them): median %d per cell (range %d-%d)\n",
            as.integer(median(n_tied)), min(n_tied), max(n_tied)))

score_set <- function(genes) {
  n_vec <- as.numeric(rownames(X) %in% genes)
  s <- vapply(seq_len(ncol(X)), function(j) {
    o <- rank_order[, j]
    wcmod.p(p_vec[o], n_vec[o], ALPHA)
  }, numeric(1))
  setNames(-s, colnames(X))
}

scores <- imap_dfr(SETS, function(genes, set_name) {
  s <- score_set(genes)
  cat(sprintf("%-22s %5d genes  score median %.2f\n", set_name, length(genes), median(s)))
  tibble(method = "miReact", target_set = set_name, cell = names(s), score = s)
})

saveRDS(scores, file.path(OUT_DIR, "scores_mireact.rds"))
cat(sprintf("Wrote scores_mireact.rds (%d target sets x %d cells)\n", length(SETS), ncol(X)))
