library(tidyverse)

BASE_DIR <- "/cfs/klemming/projects/supr/naiss2024-6-235/miRNA_footprint_new"
TOOLS_DIR <- file.path(BASE_DIR, "tools")
PROCESSED_DIR <- file.path(BASE_DIR, "datasets/HEK_SS3/processed")
OUT_DIR <- file.path(BASE_DIR, "analysis/HEK_SS3/comparisons")
MITEA_INPUT_DIR <- file.path(OUT_DIR, "mitea_input")

GFP_UPPER_PCT <- 0.99
N_TOP_TARGETS <- 200

source(file.path(BASE_DIR, "scripts/Utils.R"))  # <-- same unverified placeholder as other scripts

# --- Data --------------------------------------------------------------------
counts <- readRDS(file.path(PROCESSED_DIR, "scran_normalized_linear.rds"))
counts <- as.matrix(counts)

cell_metadata <- readRDS(file.path(PROCESSED_DIR, "cell_metadata.rds"))
GFP_counts <- cell_metadata[!is.na(cell_metadata$GFP_normalized), c("cell_id", "GFP_normalized")]
colnames(GFP_counts) <- c("Cell_ID", "eGFP")
gfp_upper_cutoff <- quantile(GFP_counts$eGFP, GFP_UPPER_PCT)

targetscan_all <- readRDS(file.path(BASE_DIR, "resources/Targets__combined_124-3p_124-3p.2_506-3p.rds"))
targetscan_sorted <- targetscan_all[order(targetscan_all$Cumulative.weighted.context...score), ]
targets_top <- intersect(targetscan_sorted[1:N_TOP_TARGETS, ]$ensembl_gene_id, rownames(counts))

# --- Our method reference ----------------------------------------------------
our_activity <- calculate_activity(counts = counts, targets = targets_top)
our_activity <- data.frame(activity = our_activity) |> rownames_to_column("Cell_ID")
our_cor_df <- merge(our_activity, GFP_counts, by = "Cell_ID")
our_cor_df <- our_cor_df[our_cor_df$eGFP < gfp_upper_cutoff, ]
our_cor_df <- our_cor_df[is.finite(our_cor_df$activity) & is.finite(our_cor_df$eGFP), ]
our_cor <- cor.test(our_cor_df$activity, our_cor_df$eGFP, method = "spearman", exact = FALSE)
cat(sprintf("Our method (scran, TargetScan): rho = %.3f, n = %d\n", our_cor$estimate, nrow(our_cor_df)))

# --- miTEA-HiRes evaluation, all four target-list variants ------------------

# --- Shared helper: read one miTEA result CSV, correlate vs GFP, sign-check -
# Factored out now that this is about to become four near-identical blocks
# (miRTarBase, TarBase, TargetScan-all, TargetScan-top200) -- same
# duplication-avoidance principle used throughout this project.
evaluate_mitea_result <- function(csv_file, score_col, label) {
  scores <- read.csv(file.path(MITEA_INPUT_DIR, csv_file))
  cor_df <- merge(scores, GFP_counts, by = "Cell_ID")
  cor_df <- cor_df[cor_df$eGFP < gfp_upper_cutoff, ]
  cor_df <- cor_df[is.finite(cor_df[[score_col]]) & is.finite(cor_df$eGFP), ]
  mitea_cor <- cor.test(cor_df[[score_col]], cor_df$eGFP, method = "spearman", exact = FALSE)
  cat(sprintf("miTEA-HiRes (scran, %s): rho = %.3f, p = %.3e, n = %d\n",
              label, mitea_cor$estimate, mitea_cor$p.value, nrow(cor_df)))

  common <- merge(our_activity, scores, by = "Cell_ID")
  sign_check <- cor.test(common$activity, common[[score_col]], method = "spearman", exact = FALSE)
  cat(sprintf("  Sanity check -- our method vs miTEA-HiRes (%s) (should be POSITIVE): rho = %.3f\n",
              label, sign_check$estimate))
  if (sign_check$estimate < 0) {
    cat(sprintf("  *** WARNING: negative sign-check for %s -- investigate before trusting this result. ***\n", label))
  }

  list(cor_df = cor_df, cor = mitea_cor, sign_check_rho = sign_check$estimate)
}

mirtarbase_res <- evaluate_mitea_result("mitea_activity_scores_mirtarbase_scran.csv", "mitea_activity", "miRTarBase")
tarbase_res <- evaluate_mitea_result("mitea_activity_scores_tarbase_scran.csv", "mitea_activity_tarbase", "TarBase")
targetscan_all_res <- evaluate_mitea_result("mitea_activity_scores_targetscan_all_scran.csv", "mitea_activity_targetscan_all", "TargetScan (all)")
targetscan_top200_res <- evaluate_mitea_result("mitea_activity_scores_targetscan_top200_scran.csv", "mitea_activity_targetscan_top200", "TargetScan (top 200)")

cat(sprintf("\n=== Comparison (scran-normalized throughout) ===\n"))
cat(sprintf("Our method, TargetScan (top 200):        rho = %.3f (n=%d)\n", our_cor$estimate, nrow(our_cor_df)))
cat(sprintf("miTEA-HiRes, miRTarBase (native):         rho = %.3f (n=%d)\n", mirtarbase_res$cor$estimate, nrow(mirtarbase_res$cor_df)))
cat(sprintf("miTEA-HiRes, TarBase:                     rho = %.3f (n=%d)\n", tarbase_res$cor$estimate, nrow(tarbase_res$cor_df)))
cat(sprintf("miTEA-HiRes, TargetScan (all candidates):  rho = %.3f (n=%d)\n", targetscan_all_res$cor$estimate, nrow(targetscan_all_res$cor_df)))
cat(sprintf("miTEA-HiRes, TargetScan (top 200):         rho = %.3f (n=%d)\n", targetscan_top200_res$cor$estimate, nrow(targetscan_top200_res$cor_df)))
cat("\nThe TargetScan (top 200) row above is the fairest direct comparison to our own method --\n",
    "same target list, same size, different algorithm. TargetScan (all) parallels the earlier\n",
    "finding that our OWN method's accuracy degraded badly on TarBase's large unranked list --\n",
    "worth checking whether miTEA-HiRes shows the same size-sensitivity or not.\n")
cat("\nFor reference, also compare against res_mireact_comparison_scran.rds (TarBase mode),\n",
    "res_bayesreact_comparison_scran.rds, and res_mireact_default_comparison_scran.rds.\n")

saveRDS(list(our = our_cor_df,
             mitea_mirtarbase = mirtarbase_res$cor_df, mitea_tarbase = tarbase_res$cor_df,
             mitea_targetscan_all = targetscan_all_res$cor_df, mitea_targetscan_top200 = targetscan_top200_res$cor_df,
             our_cor = our_cor,
             mitea_cor_mirtarbase = mirtarbase_res$cor, mitea_cor_tarbase = tarbase_res$cor,
             mitea_cor_targetscan_all = targetscan_all_res$cor, mitea_cor_targetscan_top200 = targetscan_top200_res$cor,
             sign_check_rho_mirtarbase = mirtarbase_res$sign_check_rho,
             sign_check_rho_tarbase = tarbase_res$sign_check_rho,
             sign_check_rho_targetscan_all = targetscan_all_res$sign_check_rho,
             sign_check_rho_targetscan_top200 = targetscan_top200_res$sign_check_rho,
             normalization = "scran"),
        file.path(OUT_DIR, "res_mitea_comparison_scran.rds"))
