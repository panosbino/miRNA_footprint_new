# ==============================================================================
# Step 2 (all datasets): score every target list -- benchmarked, negative
# control and random -- with our method's five variants and with miReact.
#
#   Rscript 2_score_target_lists.R <mESC|HEK|PSCSR>
#
# Our method (higher = more miRNA activity), as in "Plan for benchmarking":
#   A   raw sum                   -log2(sum of targets / its mean over cells) + 1
#   B   equal weight per gene     same, after dividing each gene by its mean
#   C   mean log ratio            -mean over targets of log2((x + 1) / (gene mean + 1))
#   C2  C minus background        background = genes in no benchmarked list of the family
#   D   mean within-cell rank     -mean over targets of the within-cell rank of the log ratios
# miReact: its TarBase-mode statistic (constant motif p = 0.01, 0/1 target
#   membership, wcmod.p() on genes ordered by expression minus gene median),
#   which scores any gene list the way miReact scores TarBase targets.
#
# Output: <out_dir>/scores/<method>.rds -- a matrix, one row per list_id, one
#   column per cell. Empty lists are NA.
# Env: N_CORES (default SLURM_CPUS_PER_TASK or 1), PSEUDOCOUNT (default 1),
#      METHODS (default all; e.g. METHODS=miReact or METHODS=Ours-A,Ours-C2 to rerun only those)
# ==============================================================================
source(file.path(Sys.getenv("MIRNA_BASE_DIR", "/cfs/klemming/projects/supr/naiss2024-6-235/miRNA_footprint_new"),
                 "scripts/benchmark_v3/benchmark_config.R"))
source(file.path(BASE_DIR, "scripts/benchmark_v3/benchmark_lib.R"))
source(file.path(MIREACT_DIR, "code/wcmod.p.R"))
ds <- dataset_from_args()
N_CORES     <- as.integer(Sys.getenv("N_CORES", Sys.getenv("SLURM_CPUS_PER_TASK", "1")))
PSEUDOCOUNT <- as.numeric(Sys.getenv("PSEUDOCOUNT", "1"))
ALL_METHODS <- c(paste0("Ours-", c("A", "B", "C", "C2", "D")), "miReact")
METHODS <- strsplit(Sys.getenv("METHODS", paste(ALL_METHODS, collapse = ",")), ",")[[1]]
if (!all(METHODS %in% ALL_METHODS)) stop("METHODS must be among: ", paste(ALL_METHODS, collapse = ", "))
score_dir <- file.path(ds$out_dir, "scores"); dir.create(score_dir, showWarnings = FALSE)

tl <- readRDS(file.path(ds$out_dir, "target_lists.rds"))
X  <- read_count_matrix(tl$expression)[tl$universe, ]
lists <- tl$lists
cat(sprintf("%s: %d lists (%s) x %d cells, %d genes; %d cores\n", ds$name, nrow(lists),
            paste(names(table(lists$kind)), table(lists$kind), sep = " ", collapse = ", "), ncol(X), nrow(X), N_CORES))

# ---- Shared transforms -------------------------------------------------------------------------
gene_mean <- rowMeans(X)
stopifnot(all(gene_mean > 0))
log_ratio <- log2((X + PSEUDOCOUNT) / (gene_mean + PSEUDOCOUNT))
log_ratio_rank <- apply(log_ratio, 2, rank, ties.method = "average")
background_mean <- map(tl$background, \(g) colMeans(log_ratio[g, , drop = FALSE]))

score_ours <- function(genes, variant, family) {
  switch(variant,
    A  = { s <- colSums(X[genes, , drop = FALSE]); -log2(s / mean(s)) + 1 },
    B  = { s <- colSums(X[genes, , drop = FALSE] / gene_mean[genes]); -log2(s / mean(s)) + 1 },
    C  = -colMeans(log_ratio[genes, , drop = FALSE]),
    C2 = -(colMeans(log_ratio[genes, , drop = FALSE]) - background_mean[[family]]),
    D  = -colMeans(log_ratio_rank[genes, , drop = FALSE]))
}

# miReact ordering: per cell, genes by (expression - gene median), decreasing
gene_median <- apply(X, 1, median)
mireact_order <- apply(X, 2, function(x) order(x - gene_median, decreasing = TRUE))
constant_p <- rep(0.01, nrow(X))
score_mireact <- function(genes) {
  is_target <- as.numeric(rownames(X) %in% genes)
  -vapply(seq_len(ncol(X)), function(j) { o <- mireact_order[, j]; wcmod.p(constant_p[o], is_target[o], 1e-10) }, numeric(1))
}

# ---- Score -----------------------------------------------------------------------------------------
empty_matrix <- function() matrix(NA_real_, nrow(lists), ncol(X), dimnames = list(lists$list_id, colnames(X)))
non_empty <- which(lists$n_genes > 0)

for (v in sub("^Ours-", "", intersect(METHODS, paste0("Ours-", c("A", "B", "C", "C2", "D"))))) {
  S <- empty_matrix()
  for (i in non_empty) S[i, ] <- score_ours(lists$genes[[i]], v, lists$family[i])
  saveRDS(S, file.path(score_dir, paste0("Ours-", v, ".rds")))
  cat(sprintf("Ours-%s done\n", v))
}

if ("miReact" %in% METHODS) {
# Free what miReact does not need: every parallel worker is a copy of this session,
# and large objects end up copied into each worker as it runs.
rm(log_ratio, log_ratio_rank, background_mean); invisible(gc())
S <- empty_matrix()
res <- parallel::mclapply(lists$genes[non_empty], score_mireact, mc.cores = N_CORES)
# A worker that dies (usually killed for using too much memory) returns an error
# instead of scores; those lists are scored again here, one by one.
failed <- which(!vapply(res, \(r) is.numeric(r) && length(r) == ncol(X), logical(1)))
if (length(failed)) {
  cat(sprintf("%d of %d lists were not returned by the parallel workers (probably out of memory);",
              length(failed), length(non_empty)),
      "scoring them again in this process. Next time, request more memory or use fewer N_CORES.\n")
  res[failed] <- lapply(lists$genes[non_empty[failed]], score_mireact)
}
for (k in seq_along(non_empty)) S[non_empty[k], ] <- res[[k]]
saveRDS(S, file.path(score_dir, "miReact.rds"))
cat(sprintf("miReact done. Scores in %s\n", score_dir))
}
