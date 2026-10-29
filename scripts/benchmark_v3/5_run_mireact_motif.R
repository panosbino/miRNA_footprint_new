# ==============================================================================
# Step 5 (all datasets): miReact motif mode for the motifs the benchmark uses.
#
#   Rscript 5_run_mireact_motif.R <mESC|HEK|PSCSR>        (N_CORES cores)
#
# miReact's addMotifExpression() scores all 16,384 7-mers and does so by
# submitting one SLURM job per 20 cells, then waiting for their results. Inside
# a pipeline job those sub-jobs are rejected silently (no account), so it waits
# forever. This script instead:
#   1. runs miReact's own addSeqs() and addMotifs(), which load its precomputed
#      per-gene motif p-values and site counts for all 7-mers,
#   2. keeps only the motifs in the dataset's motif lists (benchmarked sites,
#      control site, random 7-mers),
#   3. computes miReact's per-cell statistic for them exactly as its chunk script
#      (code/me_mrs_script.R) does -- genes ordered by expression minus gene median,
#      wcmod.p(p-values, counts, alpha = 1e-10) per motif -- in parallel over cells,
#   4. replaces +/-Inf as addMotifExpression() does and flips the sign as wrapper3.R
#      does (activity = -motif statistic).
# Each motif is computed independently, so the values equal a full miReact run for
# these motifs. Only the +/-Inf replacement value can differ (miReact uses the
# maximum over all 16,384 motifs), which never changes a cell's rank within a motif.
#
# Output: <out_dir>/mireact_motif_<dataset>.rds (motifs x cells) and
#         <out_dir>/mireact_motif_job.rds, which the summary step reads.
# ==============================================================================
source(file.path(Sys.getenv("MIRNA_BASE_DIR", "/cfs/klemming/projects/supr/naiss2024-6-235/miRNA_footprint_new"),
                 "scripts/benchmark_v3/benchmark_config.R"))
ds <- dataset_from_args()
.libPaths(c(R_LIBRARY, .libPaths()))
if (!requireNamespace("Regmex", quietly = TRUE)) stop("Regmex is not installed in ", R_LIBRARY)
N_CORES <- as.integer(Sys.getenv("N_CORES", Sys.getenv("SLURM_CPUS_PER_TASK", "1")))
ALPHA <- 1e-10                                                   # as addMotifExpression()

tl <- readRDS(file.path(ds$out_dir, "target_lists.rds"))
motifs <- sort(unique(unlist(tl$motif_lists$motifs)))
expr <- as.matrix(readRDS(tl$expression))
cat(sprintf("%s: %d motifs, %d genes x %d cells, %d cores (started %s)\n",
            ds$name, length(motifs), nrow(expr), ncol(expr), N_CORES, format(Sys.time())))

# ---- 1. miReact's own sequence and motif tables -------------------------------------------
# addMotifs() reads a global 'runparameters' and paths relative to the install directory
runparameters <- list(exp = tl$expression, motifs = 7, species = DATASETS[[ds$name]]$species, seq.type = "utr3",
                      out.meonly = TRUE, verbose = FALSE, install.dir = MIREACT_DIR, mail = NULL, tarbaserun = FALSE)
owd <- setwd(MIREACT_DIR)
source("code/addSeqs.R"); source("code/addMotifs.R"); source("code/wcmod.p.R")
sco <- addSeqs(list(exp = expr), species = runparameters$species, seq.type = runparameters$seq.type)
sco <- addMotifs(sco, motifs = runparameters$motifs)
setwd(owd)

# ---- 2. keep the benchmark's motifs ---------------------------------------------------------------
row_of <- match(motifs, names(sco$motifModels))                 # rows follow names(motifModels), as in miReact
if (anyNA(row_of)) stop("Motifs not among miReact's 7-mers: ", paste(motifs[is.na(row_of)], collapse = ", "))
if (!is.null(rownames(sco$pval.mat))) stopifnot(identical(rownames(sco$pval.mat)[row_of], motifs))
pvals  <- sco$pval.mat[row_of, , drop = FALSE]
counts <- sco$motifCounts[row_of, , drop = FALSE]
X <- sco$exp                                                     # genes with sequences, in the tables' column order
stopifnot(ncol(pvals) == nrow(X), ncol(counts) == nrow(X))
gene_median <- apply(X, 1, median)
rm(sco); invisible(gc())

# ---- 3. per-cell statistic, in parallel ---------------------------------------------------------------
score_cell <- function(j) {
  o <- order(X[, j] - gene_median, decreasing = TRUE)
  vapply(seq_along(motifs), function(m) wcmod.p(pvals[m, o], counts[m, o], ALPHA), numeric(1))
}
res <- parallel::mclapply(seq_len(ncol(X)), score_cell, mc.cores = N_CORES)
failed <- which(!vapply(res, \(r) is.numeric(r) && length(r) == length(motifs), logical(1)))
if (length(failed)) {
  cat(sprintf("%d cells not returned by the parallel workers (probably memory); computing them here.\n", length(failed)))
  res[failed] <- lapply(failed, score_cell)
}
me <- do.call(cbind, res)
dimnames(me) <- list(motifs, colnames(X))

# ---- 4. as addMotifExpression() and wrapper3.R ----------------------------------------------------------
n_inf <- sum(!is.finite(me))
me[me == Inf]  <- max(me[is.finite(me)]) + 0.01
me[me == -Inf] <- min(me[is.finite(me)]) - 0.01
me <- -me                                                         # activity = -motif statistic (wrapper3.R)

out_file <- sprintf("mireact_motif_%s.rds", ds$name)
saveRDS(me, file.path(ds$out_dir, out_file))
saveRDS(list(wd = ds$out_dir, out_file = out_file, motifs = motifs, finished = Sys.time()),
        file.path(ds$out_dir, "mireact_motif_job.rds"))
cat(sprintf("Done (%s): %d motifs x %d cells, %d infinite values replaced. Wrote %s\n",
            format(Sys.time()), nrow(me), ncol(me), n_inf, file.path(ds$out_dir, out_file)))
