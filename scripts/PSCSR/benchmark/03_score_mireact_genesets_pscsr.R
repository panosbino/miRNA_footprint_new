# ==============================================================================
# PSCSR benchmark, step 3: miReact with each target set (TarBase-mode route:
# constant p = 0.01 + 0/1 membership through wcmod.p(); verified identical to
# native addSeqs()+addMotifs() in the mESC benchmark). Parallel over sets.
# Output: analysis/PSCSR/benchmark/scores_mireact.rds
# Env: N_CORES (default SLURM_CPUS_PER_TASK or 1)
# ==============================================================================
suppressPackageStartupMessages(library(tidyverse))
source(file.path(Sys.getenv("MIRNA_BASE_DIR", "/cfs/klemming/projects/supr/naiss2024-6-235/miRNA_footprint_new"),
                 "scripts/PSCSR/benchmark/config_pscsr.R"))
source(file.path(MIREACT_DIR, "code/wcmod.p.R"))
N_CORES <- as.integer(Sys.getenv("N_CORES", Sys.getenv("SLURM_CPUS_PER_TASK", "1")))

ts_obj <- readRDS(file.path(OUT_DIR, "target_sets.rds"))
SETS <- ts_obj$sets[lengths(ts_obj$sets) > 0]
X <- readRDS(file.path(DATA_DIR, "scran_normalized_linear.rds"))[ts_obj$universe, ]
cat(sprintf("X: %d genes x %d cells; %d sets; %d cores\n", nrow(X), ncol(X), length(SETS), N_CORES))

medianExp  <- apply(X, 1, median)
rank_order <- apply(X, 2, function(x) order(x - medianExp, decreasing = TRUE))
p_vec      <- rep(0.01, nrow(X))
cat(sprintf("Genes tied at their median: median %d per cell\n", as.integer(median(colSums(sweep(X, 1, medianExp) == 0)))))

score_set <- function(genes) {
  n_vec <- as.numeric(rownames(X) %in% genes)
  -vapply(seq_len(ncol(X)), function(j) { o <- rank_order[, j]; wcmod.p(p_vec[o], n_vec[o], 1e-10) }, numeric(1))
}
res <- parallel::mclapply(SETS, score_set, mc.cores = N_CORES)
scores <- imap_dfr(res, \(s, key) tibble(method = "miReact", target_set = key, cell = colnames(X), score = s))
saveRDS(scores, file.path(OUT_DIR, "scores_mireact.rds"))
cat(sprintf("Wrote scores_mireact.rds (%d sets x %d cells)\n", length(SETS), ncol(X)))
