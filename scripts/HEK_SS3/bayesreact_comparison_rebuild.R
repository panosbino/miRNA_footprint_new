library(tidyverse)

# ---------------------------------------------------------------------------
# Practical substitute for miReact's DEFAULT motif mode (which is blocked --
# see chat discussion: required precomputed per-species k-mer probability
# files aren't published anywhere, and can't be worked around just by
# having SLURM access). bayesReact's motif_prob() implements the same
# category of computation (Markov-background-model-based k-mer enrichment)
# from scratch, with no missing-file dependency -- confirmed directly from
# the bayesReact paper's own Methods text.
#
# KNOWN FIXES REAPPLIED FROM THE ORIGINAL BUILD (not rediscovered):
#   1. process_raw_input()'s hardcoded 50%-of-larger-set gene-overlap check
#      is bypassed -- mathematically unreachable here since the reference
#      3'UTR sequence database covers fewer genes than our full panel.
#      Bypassed by calling the underlying functions it wraps directly.
#   2. build_seq_list()'s dataframe pathway does NOT deduplicate multiple
#      transcripts per gene (only its FASTA pathway does) -- replicated
#      manually (keep longest transcript per gene) before use.
#   3. bayesReact_core()'s single-motif output is a plain data.frame
#      (columns: mean, activity, post_prob, sd, "10%", "90%", n_eff, Rhat),
#      NOT the nested list-of-matrices structure the multi-motif case uses.
#
# DESIGN DECISION: scran-normalized data fed through exp_type="CPM". Checked
# against norm_scale_seq()'s actual behavior (verified earlier): for
# type="CPM" it ONLY does log2(x+1), no additional library-size rescaling --
# so it doesn't assume raw-CPM input specifically, just something already
# reasonably size-normalized, which scran-normalized data satisfies equally
# well. Reasoning, not empirically verified -- same caveat as the miReact
# rebuild; treat the sign-coherence check below as the first real test.
# ---------------------------------------------------------------------------

BASE_DIR <- "/cfs/klemming/projects/snic/naiss2024-6-235/miRNA_footprint_new"
TOOLS_DIR <- file.path(BASE_DIR, "tools")
MIREACT_DIR <- file.path(TOOLS_DIR, "miReact")
PROCESSED_DIR <- file.path(BASE_DIR, "datasets/HEK_SS3/processed")
OUT_DIR <- file.path(BASE_DIR, "analysis/HEK_SS3/comparisons")
dir.create(file.path(OUT_DIR, "bayesReact_out"), recursive = TRUE, showWarnings = FALSE)

.libPaths(c(file.path(TOOLS_DIR, "R_library"), .libPaths()))
library(bayesReact)

TARGET_MIRNA <- "hsa-miR-124-3p"
SEED_MOTIF <- "GTGCCTT"   # 7mer-m8 target site, derived + cross-validated against
                          # literature earlier in this project -- see design notes
                          # from that turn if this needs re-verifying.
GFP_UPPER_PCT <- 0.99
N_TOP_TARGETS <- 200

# TODO: same unresolved placeholder as the other rebuilt scripts.
source(file.path(BASE_DIR, "scripts/HEK_SS3/Utils.R"))  # <-- UNVERIFIED PATH

# --- Data: scran-normalized counts + GFP ------------------------------------
counts <- readRDS(file.path(PROCESSED_DIR, "scran_normalized_linear.rds"))
counts <- as.matrix(counts)
cat(sprintf("Counts: %d genes x %d cells (scran-normalized)\n", nrow(counts), ncol(counts)))

cell_metadata <- readRDS(file.path(PROCESSED_DIR, "cell_metadata.rds"))
GFP_counts <- cell_metadata[!is.na(cell_metadata$GFP_normalized), c("cell_id", "GFP_normalized")]
colnames(GFP_counts) <- c("Cell_ID", "eGFP")
gfp_upper_cutoff <- quantile(GFP_counts$eGFP, GFP_UPPER_PCT)

targetscan_all <- readRDS(file.path(BASE_DIR, "resources/Targets__combined_124-3p_124-3p.2_506-3p.rds"))
targetscan_sorted <- targetscan_all[order(targetscan_all$Cumulative.weighted.context...score), ]
targets_top <- intersect(targetscan_sorted[1:N_TOP_TARGETS, ]$ensembl_gene_id, rownames(counts))

# --- Our method reference, same data, for comparison -----------------------
our_activity <- calculate_activity(counts = counts, targets = targets_top)
our_activity <- data.frame(activity = our_activity) |> rownames_to_column("Cell_ID")
our_cor_df <- merge(our_activity, GFP_counts, by = "Cell_ID")
our_cor_df <- our_cor_df[our_cor_df$eGFP < gfp_upper_cutoff, ]
our_cor_df <- our_cor_df[is.finite(our_cor_df$activity) & is.finite(our_cor_df$eGFP), ]
our_cor <- cor.test(our_cor_df$activity, our_cor_df$eGFP, method = "spearman", exact = FALSE)
cat(sprintf("Our method (scran, TargetScan): rho = %.3f, n = %d\n", our_cor$estimate, nrow(our_cor_df)))

# --- Sequence data, deduplicated (fix #2) ------------------------------------
raw_seqs <- readRDS(file.path(MIREACT_DIR, "seqs", "hs.utr3.seqs.rds"))
seqs_deduped <- raw_seqs %>%
  group_by(gid) %>%
  filter(nchar == max(nchar)) %>%
  slice(1) %>%
  ungroup() %>%
  select(gid, sequence, nchar) %>%
  as.data.frame()
stopifnot(identical(colnames(seqs_deduped), c("gid", "sequence", "nchar")))
stopifnot(anyDuplicated(seqs_deduped$gid) == 0)
cat(sprintf("Deduplicated sequence data: %d genes (one row per gene)\n", nrow(seqs_deduped)))

# --- Bypass process_raw_input()'s 50% check (fix #1) ------------------------
exp_norm <- bayesReact::norm_scale_seq(counts, data_type = "CPM", save_rds = FALSE)

gene_set <- intersect(seqs_deduped$gid, rownames(exp_norm))
cat(sprintf("Matched genes: %d (%.1f%% of seqs_deduped, %.1f%% of counts matrix)\n",
            length(gene_set), 100 * length(gene_set) / nrow(seqs_deduped), 100 * length(gene_set) / nrow(exp_norm)))

exp_matched <- exp_norm[gene_set, , drop = FALSE]
seqs_matched <- seqs_deduped[seqs_deduped$gid %in% gene_set, ]
seqs_matched <- seqs_matched[match(gene_set, seqs_matched$gid), ]

seqlist_out <- bayesReact::build_seq_list(seqs_matched, gene_id = "gid")
motif_paths <- bayesReact::motif_prob(SEED_MOTIF, seqlist_out$seqs, seqlist_out$seqlist,
                                       paths = FALSE, cores = parallel::detectCores(),
                                       out_path = file.path(OUT_DIR, "bayesReact_out/"), include_counts = TRUE)
FC_rank_path <- bayesReact::rank_seq(exp_matched, data_type = "norm_scale_exp",
                                      path = file.path(OUT_DIR, "bayesReact_out/"))

out_paths <- list(FC_rank_path = FC_rank_path,
                   motif_probs_path = motif_paths$motif_probs_path,
                   motif_counts_path = motif_paths$motif_counts_path)

# --- Core inference -----------------------------------------------------------
bayesreact_result <- bayesReact_core(
  lst_data = list(FC_rank = out_paths$FC_rank_path,
                   motif_probs = out_paths$motif_probs_path,
                   motif_counts = out_paths$motif_counts_path),
  model = "bayesReact",
  output_type = "activity_summary",
  CI = c(0.1, 0.9)
)

str(bayesreact_result)

# --- Extraction, using the CORRECT structure (fix #3) -----------------------
bayesreact_activity <- data.frame(bayesreact_activity = bayesreact_result$activity,
                                   row.names = rownames(bayesreact_result))
ci_width <- bayesreact_result[["90%"]] - bayesreact_result[["10%"]]
cat(sprintf("bayesReact 80%% CI width: median %.3f, range [%.3f, %.3f]\n",
            median(ci_width, na.rm = TRUE), min(ci_width, na.rm = TRUE), max(ci_width, na.rm = TRUE)))
cat(sprintf("MCMC diagnostics: median n_eff = %.0f, max Rhat = %.4f\n",
            median(bayesreact_result$n_eff, na.rm = TRUE), max(bayesreact_result$Rhat, na.rm = TRUE)))

cor_df <- merge(bayesreact_activity |> rownames_to_column("Cell_ID"), GFP_counts, by = "Cell_ID")
cor_df <- cor_df[cor_df$eGFP < gfp_upper_cutoff, ]
cor_df <- cor_df[is.finite(cor_df$bayesreact_activity) & is.finite(cor_df$eGFP), ]

bayesreact_cor <- cor.test(cor_df$bayesreact_activity, cor_df$eGFP, method = "spearman", exact = FALSE)
cat(sprintf("bayesReact (scran, motif GTGCCTT): rho = %.3f, p = %.3e, n = %d\n",
            bayesreact_cor$estimate, bayesreact_cor$p.value, nrow(cor_df)))

# --- Sign-coherence check: first real test of the scran-through-CPM-pathway
common <- merge(our_activity, bayesreact_activity |> rownames_to_column("Cell_ID"), by = "Cell_ID")
sign_check <- cor.test(common$activity, common$bayesreact_activity, method = "spearman", exact = FALSE)
cat(sprintf("\nSanity check -- our method vs bayesReact (should be POSITIVE): rho = %.3f\n", sign_check$estimate))
if (sign_check$estimate < 0) {
  cat("*** WARNING: negative -- the scran-through-exp_type='CPM' reasoning may not hold. Investigate before trusting bayesreact_cor. ***\n")
}

cat(sprintf("\n=== Comparison (scran-normalized) ===\nOur method:  rho = %.3f (n=%d)\nbayesReact:  rho = %.3f (n=%d)\n",
            our_cor$estimate, nrow(our_cor_df), bayesreact_cor$estimate, nrow(cor_df)))

saveRDS(list(our = our_cor_df, bayesreact = cor_df, our_cor = our_cor, bayesreact_cor = bayesreact_cor,
             seed_motif_used = SEED_MOTIF, sign_check_rho = sign_check$estimate, normalization = "scran"),
        file.path(OUT_DIR, "res_bayesreact_comparison_scran.rds"))
