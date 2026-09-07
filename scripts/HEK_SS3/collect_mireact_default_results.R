library(tidyverse)

# ---------------------------------------------------------------------------
# Run this ONLY after confirming the launched job has actually finished
# (check `squeue -u $USER` shows nothing running, or look for the output
# file directly). Reads the tracking info saved by
# run_mireact_default_launch.R to find the timestamped job directory.
# ---------------------------------------------------------------------------

BASE_DIR <- "/cfs/klemming/projects/snic/naiss2024-6-235/miRNA_footprint_new"
PROCESSED_DIR <- file.path(BASE_DIR, "datasets/HEK_SS3/processed")
OUT_DIR <- file.path(BASE_DIR, "analysis/HEK_SS3/comparisons")

TARGET_MIRNA <- "hsa-miR-124-3p"
SEED_MOTIF <- "GTGCCTT"   # same motif derived/cross-validated earlier for bayesReact;
                          # confirmed as a valid row key below, not assumed.
GFP_UPPER_PCT <- 0.99
N_TOP_TARGETS <- 200

source(file.path(BASE_DIR, "scripts/HEK_SS3/Utils.R"))  # <-- same unverified placeholder as other scripts

tracking_path <- file.path(OUT_DIR, "mireact_default_job_tracking.rds")
if (!file.exists(tracking_path)) {
  stop(sprintf("No tracking file found at %s. Run run_mireact_default_launch.R first.", tracking_path))
}
tracking_info <- readRDS(tracking_path)
cat(sprintf("Job submitted at: %s\nExpected directory: %s\n", tracking_info$submitted_at, tracking_info$wd))

result_file <- file.path(tracking_info$wd, tracking_info$out_file)
if (!file.exists(result_file)) {
  stop(sprintf("Result file not found at %s -- the job may still be running (check `squeue -u $USER`) ",
               "or may have failed (check %s/Rscript-*.out for errors).", result_file, tracking_info$wd))
}
cat(sprintf("Found result file: %s\n", result_file))

# --- Load the motif-activity matrix -----------------------------------------
ma <- readRDS(result_file)
cat(sprintf("Motif-activity matrix: %d motifs x %d cells\n", nrow(ma), ncol(ma)))

# Sanity check on column names (the original tutorial needed a cleanup step
# here for its own data -- check ours directly rather than assume it's needed)
cat("Sample column names (should look like real cell IDs):\n")
print(head(colnames(ma)))

# --- Mandatory sanity check: is our seed motif actually a valid row key? ---
if (!(SEED_MOTIF %in% rownames(ma))) {
  stop(sprintf("'%s' not found in rownames(ma). Motif matrices are indexed by literal 7-mer ",
               "sequence -- check rownames(ma) directly (e.g. head(rownames(ma))) to see the ",
               "actual format before assuming this motif string is wrong.", SEED_MOTIF))
}
mireact_default_score <- ma[SEED_MOTIF, ]
cat(sprintf("*** Confirmed '%s' present as a row in ma. Score range: [%.3f, %.3f] ***\n",
            SEED_MOTIF, min(mireact_default_score), max(mireact_default_score)))

mireact_df <- data.frame(mireact_default_activity = mireact_default_score) |> rownames_to_column("Cell_ID")

# --- GFP + our own method's reference, same convention as other rebuilt scripts
cell_metadata <- readRDS(file.path(PROCESSED_DIR, "cell_metadata.rds"))
GFP_counts <- cell_metadata[!is.na(cell_metadata$GFP_normalized), c("cell_id", "GFP_normalized")]
colnames(GFP_counts) <- c("Cell_ID", "eGFP")
gfp_upper_cutoff <- quantile(GFP_counts$eGFP, GFP_UPPER_PCT)

counts <- readRDS(file.path(PROCESSED_DIR, "scran_normalized_linear.rds"))
counts <- as.matrix(counts)
targetscan_all <- readRDS(file.path(BASE_DIR, "resources/Targets__combined_124-3p_124-3p.2_506-3p.rds"))
targetscan_sorted <- targetscan_all[order(targetscan_all$Cumulative.weighted.context...score), ]
targets_top <- intersect(targetscan_sorted[1:N_TOP_TARGETS, ]$ensembl_gene_id, rownames(counts))

our_activity <- calculate_activity(counts = counts, targets = targets_top)
our_activity <- data.frame(activity = our_activity) |> rownames_to_column("Cell_ID")
our_cor_df <- merge(our_activity, GFP_counts, by = "Cell_ID")
our_cor_df <- our_cor_df[our_cor_df$eGFP < gfp_upper_cutoff, ]
our_cor_df <- our_cor_df[is.finite(our_cor_df$activity) & is.finite(our_cor_df$eGFP), ]
our_cor <- cor.test(our_cor_df$activity, our_cor_df$eGFP, method = "spearman", exact = FALSE)

# --- miReact default-mode result --------------------------------------------
cor_df <- merge(mireact_df, GFP_counts, by = "Cell_ID")
cor_df <- cor_df[cor_df$eGFP < gfp_upper_cutoff, ]
cor_df <- cor_df[is.finite(cor_df$mireact_default_activity) & is.finite(cor_df$eGFP), ]
mireact_default_cor <- cor.test(cor_df$mireact_default_activity, cor_df$eGFP, method = "spearman", exact = FALSE)
cat(sprintf("miReact (default motif-search mode, scran): rho = %.3f, p = %.3e, n = %d\n",
            mireact_default_cor$estimate, mireact_default_cor$p.value, nrow(cor_df)))

# --- Sign-coherence check ----------------------------------------------------
common <- merge(our_activity, mireact_df, by = "Cell_ID")
sign_check <- cor.test(common$activity, common$mireact_default_activity, method = "spearman", exact = FALSE)
cat(sprintf("\nSanity check -- our method vs miReact default-mode (should be POSITIVE): rho = %.3f\n",
            sign_check$estimate))

cat(sprintf("\n=== Comparison (scran-normalized) ===\n"))
cat(sprintf("Our method, TargetScan:        rho = %.3f (n=%d)\n", our_cor$estimate, nrow(our_cor_df)))
cat(sprintf("miReact, default motif-search: rho = %.3f (n=%d)\n", mireact_default_cor$estimate, nrow(cor_df)))
cat("\nFor reference, also compare against the earlier TarBase-mode and bayesReact motif-mode results\n",
    "(res_mireact_comparison_scran.rds, res_bayesreact_comparison_scran.rds) -- three genuinely\n",
    "different methods/target-definitions now available for the same miRNA on the same data.\n")

saveRDS(list(our = our_cor_df, mireact_default = cor_df, our_cor = our_cor,
             mireact_default_cor = mireact_default_cor, seed_motif_used = SEED_MOTIF,
             sign_check_rho = sign_check$estimate, normalization = "scran"),
        file.path(OUT_DIR, "res_mireact_default_comparison_scran.rds"))
