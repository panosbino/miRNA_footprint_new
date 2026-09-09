library(tidyverse)

BASE_DIR <- "/cfs/klemming/projects/supr/naiss2024-6-235/miRNA_footprint_new"
TOOLS_DIR <- file.path(BASE_DIR, "tools")
PROCESSED_DIR <- file.path(BASE_DIR, "datasets/HEK_SS3/processed")
OUT_DIR <- file.path(BASE_DIR, "analysis/HEK_SS3/comparisons")
MITEA_INPUT_DIR <- file.path(OUT_DIR, "mitea_input")

GFP_UPPER_PCT <- 0.99
N_TOP_TARGETS <- 200

source(file.path(BASE_DIR, "scripts/HEK_SS3/Utils.R"))  # <-- same unverified placeholder as other scripts

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

# --- miTEA-HiRes, native miRTarBase -------------------------------------------
mitea_scores <- read.csv(file.path(MITEA_INPUT_DIR, "mitea_activity_scores_mirtarbase_scran.csv"))
cor_df <- merge(mitea_scores, GFP_counts, by = "Cell_ID")
cor_df <- cor_df[cor_df$eGFP < gfp_upper_cutoff, ]
cor_df <- cor_df[is.finite(cor_df$mitea_activity) & is.finite(cor_df$eGFP), ]
mitea_cor <- cor.test(cor_df$mitea_activity, cor_df$eGFP, method = "spearman", exact = FALSE)
cat(sprintf("miTEA-HiRes (scran, miRTarBase): rho = %.3f, p = %.3e, n = %d\n",
            mitea_cor$estimate, mitea_cor$p.value, nrow(cor_df)))

# --- miTEA-HiRes, TarBase-matched --------------------------------------------
mitea_scores_tb <- read.csv(file.path(MITEA_INPUT_DIR, "mitea_activity_scores_tarbase_scran.csv"))
cor_df_tb <- merge(mitea_scores_tb, GFP_counts, by = "Cell_ID")
cor_df_tb <- cor_df_tb[cor_df_tb$eGFP < gfp_upper_cutoff, ]
cor_df_tb <- cor_df_tb[is.finite(cor_df_tb$mitea_activity_tarbase) & is.finite(cor_df_tb$eGFP), ]
mitea_cor_tb <- cor.test(cor_df_tb$mitea_activity_tarbase, cor_df_tb$eGFP, method = "spearman", exact = FALSE)
cat(sprintf("miTEA-HiRes (scran, TarBase): rho = %.3f, p = %.3e, n = %d\n",
            mitea_cor_tb$estimate, mitea_cor_tb$p.value, nrow(cor_df_tb)))

# --- Sign-coherence checks (first real test of scran-through-normalize_counts()) --
common <- merge(our_activity, mitea_scores, by = "Cell_ID")
sign_check <- cor.test(common$activity, common$mitea_activity, method = "spearman", exact = FALSE)
cat(sprintf("\nSanity check -- our method vs miTEA-HiRes (miRTarBase) (should be POSITIVE): rho = %.3f\n",
            sign_check$estimate))

common_tb <- merge(our_activity, mitea_scores_tb, by = "Cell_ID")
sign_check_tb <- cor.test(common_tb$activity, common_tb$mitea_activity_tarbase, method = "spearman", exact = FALSE)
cat(sprintf("Sanity check -- our method vs miTEA-HiRes (TarBase) (should be POSITIVE): rho = %.3f\n",
            sign_check_tb$estimate))

if (sign_check$estimate < 0 || sign_check_tb$estimate < 0) {
  cat("*** WARNING: negative correlation found -- the scran-through-normalize_counts() reasoning\n",
      "may not hold. Investigate before trusting mitea_cor / mitea_cor_tb as-is. ***\n")
}

cat(sprintf("\n=== Comparison (scran-normalized throughout) ===\n"))
cat(sprintf("Our method, TargetScan:      rho = %.3f (n=%d)\n", our_cor$estimate, nrow(our_cor_df)))
cat(sprintf("miTEA-HiRes, miRTarBase:     rho = %.3f (n=%d)\n", mitea_cor$estimate, nrow(cor_df)))
cat(sprintf("miTEA-HiRes, TarBase:        rho = %.3f (n=%d)\n", mitea_cor_tb$estimate, nrow(cor_df_tb)))
cat("\nFor reference, compare against res_mireact_comparison_scran.rds (TarBase mode),\n",
    "res_bayesreact_comparison_scran.rds, and res_mireact_default_comparison_scran.rds --\n",
    "four methods now available on the same scran-normalized data.\n")

saveRDS(list(our = our_cor_df, mitea_mirtarbase = cor_df, mitea_tarbase = cor_df_tb,
             our_cor = our_cor, mitea_cor_mirtarbase = mitea_cor, mitea_cor_tarbase = mitea_cor_tb,
             sign_check_rho = sign_check$estimate, sign_check_tb_rho = sign_check_tb$estimate,
             normalization = "scran"),
        file.path(OUT_DIR, "res_mitea_comparison_scran.rds"))
