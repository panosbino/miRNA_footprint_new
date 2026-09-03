library(tidyverse)

# ---------------------------------------------------------------------------
# REBUILT from scratch under the new miRNA_footprint_new structure. Run
# setup_tools.sh once first (clones miReact, installs Regmex).
#
# DESIGN DECISION, stated explicitly rather than left implicit: both our
# method AND miReact use scran-normalized counts here. miReact's scoring is
# rank-based (ranks genes within each cell relative to that cell's median
# expression), which is largely invariant to WHICH monotonic per-cell
# normalization was used -- so feeding it scran-normalized data instead of
# its usual convention is a defensible, low-risk choice, unlike bayesReact
# (log2(x+1)-specific) or miTEA-HiRes (does its own internal normalization
# regardless of input). Those two keep their native processing until
# separately verified -- do not assume this same reasoning extends to them.
#
# KNOWN BUGS IN MIREACT'S OWN SOURCE, reapplied from prior discovery (all
# reconfirmed present in a fresh clone before writing this):
#   - addMotifs() reads runparameters$tarbaserun AND runparameters$species
#     from the GLOBAL environment, not sco$runparameters -- set both
#     explicitly.
#   - `tar <- readRDS(tar, file="./data/tarbase.rds")` passes an undefined
#     `tar` positionally alongside a named file= -- appears harmless (lazy
#     evaluation likely means the positional arg, bound to refhook, is
#     never forced) but flagged in case behavior differs on this R version.
# ---------------------------------------------------------------------------

BASE_DIR <- "/cfs/klemming/projects/snic/naiss2024-6-235/miRNA_footprint_new"
TOOLS_DIR <- file.path(BASE_DIR, "tools")
MIREACT_DIR <- file.path(TOOLS_DIR, "miReact")
PROCESSED_DIR <- file.path(BASE_DIR, "datasets/HEK_SS3/processed")
OUT_DIR <- file.path(BASE_DIR, "analysis/HEK_SS3/comparisons")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

.libPaths(c(file.path(TOOLS_DIR, "R_library"), .libPaths()))

TARGET_MIRNA <- "hsa-miR-124-3p"
GFP_UPPER_PCT <- 0.99
N_TOP_TARGETS <- 200

# TODO: PLACEHOLDER, same unresolved item as depth_subsampling_scran_new.R --
# confirm the real Utils.R path before running.
source(file.path(BASE_DIR, "scripts/HEK_SS3/Utils.R"))  # <-- UNVERIFIED PATH

# --- Data: scran-normalized counts (precomputed, full dataset) -------------
counts <- readRDS(file.path(PROCESSED_DIR, "scran_normalized_linear.rds"))
counts <- as.matrix(counts)
cat(sprintf("Counts: %d genes x %d cells (scran-normalized, protein-coding-restricted)\n",
            nrow(counts), ncol(counts)))

cell_metadata <- readRDS(file.path(PROCESSED_DIR, "cell_metadata.rds"))
stopifnot(all(c("cell_id", "GFP_normalized") %in% colnames(cell_metadata)))
GFP_counts <- cell_metadata[!is.na(cell_metadata$GFP_normalized), c("cell_id", "GFP_normalized")]
colnames(GFP_counts) <- c("Cell_ID", "eGFP")
gfp_upper_cutoff <- quantile(GFP_counts$eGFP, GFP_UPPER_PCT)
cat(sprintf("GFP: %d cells with non-NA GFP_normalized; fixed QC cutoff (%.0fth pct) = %.3f\n",
            nrow(GFP_counts), GFP_UPPER_PCT * 100, gfp_upper_cutoff))

targetscan_all <- readRDS(file.path(BASE_DIR, "resources/Targets__combined_124-3p_124-3p.2_506-3p.rds"))
targetscan_sorted <- targetscan_all[order(targetscan_all$Cumulative.weighted.context...score), ]
targets_top <- intersect(targetscan_sorted[1:N_TOP_TARGETS, ]$ensembl_gene_id, rownames(counts))
cat(sprintf("Of top %d TargetScan targets, %d present in scran-normalized matrix (%d missing)\n",
            N_TOP_TARGETS, length(targets_top), N_TOP_TARGETS - length(targets_top)))

# --- Our method (aggregative), scran-normalized, TargetScan targets --------
our_activity <- calculate_activity(counts = counts, targets = targets_top)
our_activity <- data.frame(activity = our_activity) |> rownames_to_column("Cell_ID")
our_cor_df <- merge(our_activity, GFP_counts, by = "Cell_ID")
our_cor_df <- our_cor_df[our_cor_df$eGFP < gfp_upper_cutoff, ]
our_cor_df <- our_cor_df[is.finite(our_cor_df$activity) & is.finite(our_cor_df$eGFP), ]
our_cor <- cor.test(our_cor_df$activity, our_cor_df$eGFP, method = "spearman", exact = FALSE)
cat(sprintf("Our method (scran, TargetScan): rho = %.3f, p = %.3e, n = %d\n",
            our_cor$estimate, our_cor$p.value, nrow(our_cor_df)))

# --- TarBase targets, for our own method too (matches TarBase-controlled --
# --- comparison design used earlier in this project) ------------------------
tar <- readRDS(file.path(MIREACT_DIR, "data", "tarbase.rds"))
tar_124 <- tar[tar$mirna == TARGET_MIRNA & tar$species == "Homo sapiens" &
                 tar$up_down == "DOWN" & !is.na(tar$up_down), ]
tarbase_symbols <- unique(sub("\\(hsa\\)$", "", tar_124$geneName))
cat(sprintf("TarBase %s DOWN targets: %d unique symbols\n", TARGET_MIRNA, length(tarbase_symbols)))

# --- miReact setup, reapplying known fixes ----------------------------------
setwd(MIREACT_DIR)
source("code/addSeqs.R")
source("code/addMotifs.R")
source("code/wcmod.p.R")

sco <- list(exp = counts)   # scran-normalized, per design decision above
sco <- addSeqs(sco, species = "hs", seq.type = "utr3")

runparameters <- list(tarbaserun = TRUE, species = "hs")  # global-env fix, reapplied
sco <- addMotifs(sco, motifs = 7)

mirna_idx <- which(rownames(sco$motifCounts) == TARGET_MIRNA)
if (length(mirna_idx) == 0) stop(sprintf("'%s' not found in sco$motifCounts rownames.", TARGET_MIRNA))
n_nonzero_targets <- sum(sco$motifCounts[mirna_idx, ] > 0)
cat(sprintf("*** SANITY CHECK: %d of %d genes matched as %s targets (TarBase, via miReact) ***\n",
            n_nonzero_targets, ncol(sco$motifCounts), TARGET_MIRNA))
if (n_nonzero_targets == 0) stop("Zero target matches -- do not proceed. Check symbol mapping.")

medianExp <- apply(sco$exp, 1, median)
rank_order <- apply(sco$exp, 2, function(x) order(x - medianExp, decreasing = TRUE))
alpha <- 1e-10
pval_row <- sco$pval.mat[mirna_idx, ]
counts_row <- sco$motifCounts[mirna_idx, ]

cat("Scoring", ncol(sco$exp), "cells for", TARGET_MIRNA, "...\n")
mireact_score <- sapply(seq_len(ncol(sco$exp)), function(s) {
  wcmod.p(pval_row[rank_order[, s]], counts_row[rank_order[, s]], alpha)
})
names(mireact_score) <- colnames(sco$exp)
mireact_score <- mireact_score * -1   # sign convention, per wrapper3.R

mireact_df <- data.frame(mireact_activity = mireact_score) |> rownames_to_column("Cell_ID")

# --- Our method, TarBase targets (algorithm-only comparison vs miReact) ----
setwd(OUT_DIR)  # back out of the miReact repo dir before touching our own paths
gsym_to_gid_map <- sco$seqs[, c("gid", "gsym")]
gsym_to_gid_map <- gsym_to_gid_map[!duplicated(gsym_to_gid_map$gsym), ]
tarbase_gids <- intersect(gsym_to_gid_map$gid[match(tarbase_symbols, gsym_to_gid_map$gsym)], rownames(counts))
tarbase_gids <- tarbase_gids[!is.na(tarbase_gids)]
cat(sprintf("Of %d TarBase symbols, %d mapped to a gid present in the scran-normalized matrix\n",
            length(tarbase_symbols), length(tarbase_gids)))

our_tarbase_activity <- calculate_activity(counts = counts, targets = tarbase_gids)
our_tarbase_activity <- data.frame(activity = our_tarbase_activity) |> rownames_to_column("Cell_ID")
our_tarbase_cor_df <- merge(our_tarbase_activity, GFP_counts, by = "Cell_ID")
our_tarbase_cor_df <- our_tarbase_cor_df[our_tarbase_cor_df$eGFP < gfp_upper_cutoff, ]
our_tarbase_cor_df <- our_tarbase_cor_df[is.finite(our_tarbase_cor_df$activity) & is.finite(our_tarbase_cor_df$eGFP), ]
our_tarbase_cor <- cor.test(our_tarbase_cor_df$activity, our_tarbase_cor_df$eGFP, method = "spearman", exact = FALSE)
cat(sprintf("Our method (scran, TarBase): rho = %.3f, p = %.3e, n = %d\n",
            our_tarbase_cor$estimate, our_tarbase_cor$p.value, nrow(our_tarbase_cor_df)))

# --- miReact result -----------------------------------------------------------
cor_df <- merge(mireact_df, GFP_counts, by = "Cell_ID")
cor_df <- cor_df[cor_df$eGFP < gfp_upper_cutoff, ]
cor_df <- cor_df[is.finite(cor_df$mireact_activity) & is.finite(cor_df$eGFP), ]
mireact_cor <- cor.test(cor_df$mireact_activity, cor_df$eGFP, method = "spearman", exact = FALSE)
cat(sprintf("miReact (scran, TarBase): rho = %.3f, p = %.3e, n = %d\n",
            mireact_cor$estimate, mireact_cor$p.value, nrow(cor_df)))

# --- Sign-coherence check (same principle as before: agreement on real ------
# --- targets doesn't guarantee agreement in a null/random regime, but this --
# --- isn't a random-gene test -- just a basic sanity check on real targets) -
common <- merge(our_tarbase_activity, mireact_df, by = "Cell_ID")
sign_check <- cor.test(common$activity, common$mireact_activity, method = "spearman", exact = FALSE)
cat(sprintf("\nSanity check -- our method vs miReact, both on TarBase targets (should be POSITIVE): rho = %.3f\n",
            sign_check$estimate))

cat(sprintf("\n=== Comparison (all on scran-normalized data) ===\n"))
cat(sprintf("Our method, TargetScan: rho = %.3f (n=%d)\n", our_cor$estimate, nrow(our_cor_df)))
cat(sprintf("Our method, TarBase:    rho = %.3f (n=%d)\n", our_tarbase_cor$estimate, nrow(our_tarbase_cor_df)))
cat(sprintf("miReact,    TarBase:    rho = %.3f (n=%d)\n", mireact_cor$estimate, nrow(cor_df)))

saveRDS(list(our_targetscan = our_cor_df, our_tarbase = our_tarbase_cor_df, mireact = cor_df,
             our_targetscan_cor = our_cor, our_tarbase_cor = our_tarbase_cor, mireact_cor = mireact_cor,
             n_targets_matched_mireact = n_nonzero_targets, sign_check_rho = sign_check$estimate,
             normalization = "scran"),
        file.path(OUT_DIR, "res_mireact_comparison_scran.rds"))
