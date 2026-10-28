# ==============================================================================
# Step 4 (all datasets): bayesReact for every motif in the dataset's motif lists
# (benchmarked seed sites, the control site, random 7-mers).
#
#   Rscript 4_run_bayesreact.R <mESC|HEK|PSCSR>   # all motifs, one after another (as run_benchmark_v3.sbatch does)
# With SLURM_ARRAY_TASK_ID set, only motif number <task id> runs (for splitting
# the work over a SLURM array if a single job is too slow).
#
# Same calls and fixes as the earlier bayesReact scripts: longest 3'UTR per
# gene, data_type = "CPM" (log2(x + 1)), process_raw_input() bypassed. Each
# motif works in its own folder, so array tasks never share intermediate files.
# Output: <out_dir>/bayesreact/<motif>.rds  list(activity, diagnostics(cell, n_eff, Rhat, ci_width))
# ==============================================================================
source(file.path(Sys.getenv("MIRNA_BASE_DIR", "/cfs/klemming/projects/supr/naiss2024-6-235/miRNA_footprint_new"),
                 "scripts/benchmark_v3/benchmark_config.R"))
source(file.path(BASE_DIR, "scripts/benchmark_v3/benchmark_lib.R"))
.libPaths(c(R_LIBRARY, .libPaths()))
suppressPackageStartupMessages(library(bayesReact))
ds <- dataset_from_args()
N_CORES <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", parallel::detectCores()))
out <- file.path(ds$out_dir, "bayesreact"); dir.create(out, showWarnings = FALSE)

tl <- readRDS(file.path(ds$out_dir, "target_lists.rds"))
motifs <- sort(unique(unlist(tl$motif_lists$motifs)))
task <- Sys.getenv("SLURM_ARRAY_TASK_ID", "")
cat(sprintf("%s: %d motifs to run (submit with --array=1-%d)\n", ds$name, length(motifs), length(motifs)))
todo <- if (nzchar(task)) motifs[as.integer(task)] else motifs

expr <- read_count_matrix(tl$expression)
seqs <- load_utr_table(DATASETS[[ds$name]]$species) |>
  group_by(gid) |> filter(nchar == max(nchar)) |> slice(1) |> ungroup() |>
  dplyr::select(gid, sequence, nchar) |> as.data.frame()
expr_norm <- bayesReact::norm_scale_seq(expr, data_type = "CPM", save_rds = FALSE)
genes <- intersect(seqs$gid, rownames(expr_norm))
seq_list <- bayesReact::build_seq_list(seqs[match(genes, seqs$gid), ], gene_id = "gid")
expr_norm <- expr_norm[genes, , drop = FALSE]

for (motif in todo) {
  work <- file.path(out, "work", motif); dir.create(work, recursive = TRUE, showWarnings = FALSE)
  cat(sprintf("\n=== bayesReact: motif %s ===\n", motif))
  ranks <- bayesReact::rank_seq(expr_norm, data_type = "norm_scale_exp", path = paste0(work, "/"))
  probs <- bayesReact::motif_prob(motif, seq_list$seqs, seq_list$seqlist, paths = FALSE, cores = N_CORES,
                                  out_path = paste0(work, "/"), include_counts = TRUE)
  res <- bayesReact_core(lst_data = list(FC_rank = ranks, motif_probs = probs$motif_probs_path, motif_counts = probs$motif_counts_path),
                         model = "bayesReact", output_type = "activity_summary", CI = c(0.1, 0.9))
  diagnostics <- data.frame(cell = rownames(res), n_eff = res$n_eff, Rhat = res$Rhat, ci_width = res[["90%"]] - res[["10%"]])
  cat(sprintf("MCMC: median n_eff %.0f, max Rhat %.4f, cells with Rhat > 1.01: %d\n",
              median(diagnostics$n_eff, na.rm = TRUE), max(diagnostics$Rhat, na.rm = TRUE), sum(diagnostics$Rhat > 1.01, na.rm = TRUE)))
  saveRDS(list(activity = setNames(res$activity, rownames(res)), diagnostics = diagnostics), file.path(out, paste0(motif, ".rds")))
}
