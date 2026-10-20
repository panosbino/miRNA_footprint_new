# ==============================================================================
# PSCSR benchmark, step 2: our method, variants A/B/C/C2/D x every target set.
# Same definitions as the HEK/mESC benchmarks (tarbase_diag_common.R). C2's
# background is unit-specific: genes in U in no target set of that unit.
# Output: analysis/PSCSR/benchmark/scores_ours.rds (method, target_set, cell, score)
# ==============================================================================
suppressPackageStartupMessages(library(tidyverse))
source(file.path(Sys.getenv("MIRNA_BASE_DIR", "/cfs/klemming/projects/supr/naiss2024-6-235/miRNA_footprint_new"),
                 "scripts/PSCSR/benchmark/config_pscsr.R"))
source(file.path(BASE_DIR, "scripts/Utils.R"))   # calculate_activity()
PSEUDOCOUNT <- as.numeric(Sys.getenv("PSEUDOCOUNT", "1"))

ts_obj <- readRDS(file.path(OUT_DIR, "target_sets.rds"))
SETS <- ts_obj$sets[lengths(ts_obj$sets) > 0]
X <- readRDS(file.path(DATA_DIR, "scran_normalized_linear.rds"))[ts_obj$universe, ]
cat(sprintf("X: %d genes x %d cells; %d non-empty target sets\n", nrow(X), ncol(X), length(SETS)))

gene_mean <- rowMeans(X); stopifnot(all(gene_mean > 0))
LOGRATIO  <- log2((X + PSEUDOCOUNT) / (gene_mean + PSEUDOCOUNT))
RANKS     <- apply(LOGRATIO, 2, rank, ties.method = "average")
bg_mean   <- map(ts_obj$background, \(g) colMeans(LOGRATIO[g, , drop = FALSE]))

score_variant <- function(genes, variant, unit) {
  s <- switch(variant,
    A  = calculate_activity(counts = X, targets = genes),
    B  = { v <- colSums(X[genes, , drop = FALSE] / gene_mean[genes]); -log2(v / mean(v)) + 1 },
    C  = -colMeans(LOGRATIO[genes, , drop = FALSE]),
    C2 = -(colMeans(LOGRATIO[genes, , drop = FALSE]) - bg_mean[[unit]]),
    D  = -colMeans(RANKS[genes, , drop = FALSE]))
  setNames(as.numeric(s), colnames(X))
}
VARIANTS <- c("A", "B", "C", "C2", "D")
scores <- map_dfr(VARIANTS, function(v) imap_dfr(SETS, function(genes, key) {
  s <- score_variant(genes, v, split_key(key)$unit)
  tibble(method = paste0("Ours-", v), target_set = key, cell = names(s), score = s)
}))
n_bad <- sum(!is.finite(scores$score))
if (n_bad) cat(sprintf("%d non-finite scores (+Inf = no target expression in a cell; kept)\n", n_bad))
saveRDS(scores, file.path(OUT_DIR, "scores_ours.rds"))
cat(sprintf("Wrote scores_ours.rds: %d variants x %d sets x %d cells\n", length(VARIANTS), length(SETS), ncol(X)))
