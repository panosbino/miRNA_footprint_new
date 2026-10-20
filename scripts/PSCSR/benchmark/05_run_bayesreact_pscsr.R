# ==============================================================================
# PSCSR benchmark, step 5: bayesReact, one run per unit motif (no target list).
# Same calls and fixes as mESC_bayesreact_separation.R / bayesreact_comparison_rebuild.R
# (longest 3'UTR per gene, process_raw_input() bypassed, single-motif output
# read as a data.frame), on the prepared scran matrix through data_type = "CPM"
# (= log2(x + 1)), as for HEK.
#
# One unit per call: UNIT_INDEX (1-based row of units.tsv) or SLURM_ARRAY_TASK_ID;
# if neither is set, all units run in sequence. Each unit works in its own
# folder, so array tasks cannot overwrite each other's intermediate files.
# Output: analysis/PSCSR/benchmark/bayesreact/<unit>.rds
#   list(activity = named vector, diagnostics = data.frame(cell, n_eff, Rhat, ci_width), motif)
# ==============================================================================
suppressPackageStartupMessages(library(tidyverse))
source(file.path(Sys.getenv("MIRNA_BASE_DIR", "/cfs/klemming/projects/supr/naiss2024-6-235/miRNA_footprint_new"),
                 "scripts/PSCSR/benchmark/config_pscsr.R"))
.libPaths(c(file.path(BASE_DIR, "tools/R_library"), .libPaths()))
library(bayesReact)
N_CORES <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", parallel::detectCores()))
BR_DIR  <- file.path(OUT_DIR, "bayesreact")
dir.create(BR_DIR, recursive = TRUE, showWarnings = FALSE)

units <- read_tsv(file.path(DATA_DIR, "units.tsv"), show_col_types = FALSE)
idx <- Sys.getenv("UNIT_INDEX", Sys.getenv("SLURM_ARRAY_TASK_ID", ""))
todo <- if (nzchar(idx)) as.integer(idx) else seq_len(nrow(units))
stopifnot(all(todo >= 1 & todo <= nrow(units)))

counts <- readRDS(file.path(DATA_DIR, "scran_normalized_linear.rds"))
seqs <- readRDS(file.path(MIREACT_DIR, "seqs/hs.utr3.seqs.rds")) |>
  group_by(gid) |> filter(nchar == max(nchar)) |> slice(1) |> ungroup() |>
  dplyr::select(gid, sequence, nchar) |> as.data.frame()
exp_norm <- bayesReact::norm_scale_seq(counts, data_type = "CPM", save_rds = FALSE)
gene_set <- intersect(seqs$gid, rownames(exp_norm))
exp_matched <- exp_norm[gene_set, , drop = FALSE]
seqs_matched <- seqs[match(gene_set, seqs$gid), ]
seqlist_out <- bayesReact::build_seq_list(seqs_matched, gene_id = "gid")
cat(sprintf("%d genes x %d cells; units to run: %s\n", length(gene_set), ncol(exp_matched), paste(units$unit[todo], collapse = ", ")))

for (i in todo) {
  u <- units$unit[i]; motif <- units$motif[i]
  work <- file.path(BR_DIR, "work", u); dir.create(work, recursive = TRUE, showWarnings = FALSE)
  cat(sprintf("\n=== bayesReact: %s (motif %s) ===\n", u, motif))
  FC_rank_path <- bayesReact::rank_seq(exp_matched, data_type = "norm_scale_exp", path = paste0(work, "/"))
  motif_paths <- bayesReact::motif_prob(motif, seqlist_out$seqs, seqlist_out$seqlist, paths = FALSE, cores = N_CORES,
                                        out_path = paste0(work, "/"), include_counts = TRUE)
  res <- bayesReact_core(lst_data = list(FC_rank = FC_rank_path, motif_probs = motif_paths$motif_probs_path,
                                         motif_counts = motif_paths$motif_counts_path),
                         model = "bayesReact", output_type = "activity_summary", CI = c(0.1, 0.9))
  diag <- data.frame(cell = rownames(res), n_eff = res$n_eff, Rhat = res$Rhat, ci_width = res[["90%"]] - res[["10%"]])
  cat(sprintf("MCMC: median n_eff %.0f, max Rhat %.4f, cells with Rhat > 1.01: %d\n",
              median(diag$n_eff, na.rm = TRUE), max(diag$Rhat, na.rm = TRUE), sum(diag$Rhat > 1.01, na.rm = TRUE)))
  saveRDS(list(activity = setNames(res$activity, rownames(res)), diagnostics = diag, motif = motif),
          file.path(BR_DIR, paste0(u, ".rds")))
}
