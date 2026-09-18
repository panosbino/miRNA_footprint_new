library(tidyverse)
library(Rcpp)
library(RcppArmadillo)
library(igraph)

# ---------------------------------------------------------------------------
# Data-driven target filtering: instead of relying only on TargetScan's own
# context-score ranking, use the CO-EXPRESSION structure observed directly
# in our data. Rationale: true, functionally active miR-124 targets should
# be coordinately repressed together (driven by the same underlying
# per-cell miR-124 activity), so they should correlate with EACH OTHER
# across cells -- whereas weak/false-positive predicted targets have no
# reason to show this structure. Build a graph from strong positive
# pairwise correlations among candidate targets, then keep only genes that
# are well-connected within that network (>= 3 correlating neighbors).
#
# DESIGN DECISIONS, stated explicitly:
# 1. Starting pool = FULL TargetScan candidate list, not top-200. The whole
#    point is to let co-expression do the filtering job -- pre-filtering
#    with TargetScan's own score first would partly defeat that.
# 2. Input data = scran-normalized (matches the recent focus of this
#    project). Restricted to protein-coding, zero-filtered genes, same as
#    every other scran-based script here.
# 3. Correlation threshold: positive, > 0.4 (per explicit instruction).
# 4. Degree filter: >= 3 correlating neighbors (per explicit instruction).
# ---------------------------------------------------------------------------

BASE_DIR <- "/cfs/klemming/projects/supr/naiss2024-6-235/miRNA_footprint_new"
PROCESSED_DIR <- file.path(BASE_DIR, "datasets/HEK_SS3/processed")
OUT_DIR <- file.path(BASE_DIR, "analysis/HEK_SS3/comparisons")

GFP_UPPER_PCT <- 0.99
CORR_THRESHOLD <- 0.4
MIN_DEGREE <- 3

source(file.path(BASE_DIR, "scripts/HEK_SS3/Utils.R"))  # <-- same unverified placeholder as other rebuilt scripts

# --- Compile the C++ correlation function -----------------------------------
sourceCpp(file.path(BASE_DIR, "scripts/HEK_SS3/fullSpearman.cpp"))  # adjust path to wherever you place the .cpp file

# --- MANDATORY VALIDATION: does this match R's own Spearman convention? ----
# The C++ rank function explicitly does NOT average-rank ties (confirmed by
# reading the source). Checking directly, on our own real data, whether
# this produces a meaningfully different answer than R's own cor() --
# rather than assuming either way.
cat("Validation check will run after data loading (see below) -- comparing\n",
    "spearman_full_cpp() against R's own cor(method='spearman') on a sample\n",
    "of real target gene pairs from our actual data.\n")

# --- Data --------------------------------------------------------------------
counts <- readRDS(file.path(PROCESSED_DIR, "scran_normalized_linear.rds"))
counts <- as.matrix(counts)
cat(sprintf("Counts: %d genes x %d cells (scran-normalized)\n", nrow(counts), ncol(counts)))

cell_metadata <- readRDS(file.path(PROCESSED_DIR, "cell_metadata.rds"))
GFP_counts <- cell_metadata[!is.na(cell_metadata$GFP_normalized), c("cell_id", "GFP_normalized")]
colnames(GFP_counts) <- c("Cell_ID", "eGFP")
gfp_upper_cutoff <- quantile(GFP_counts$eGFP, GFP_UPPER_PCT)

targetscan_all <- readRDS(file.path(BASE_DIR, "resources/Targets__combined_124-3p_124-3p.2_506-3p.rds"))
targetscan_ids <- intersect(targetscan_all$ensembl_gene_id, rownames(counts))
cat(sprintf("TargetScan candidates: %d total, %d present in scran-normalized matrix\n",
            nrow(targetscan_all), length(targetscan_ids)))

# --- Restrict to candidate targets, drop zero-variance genes ---------------
# Spearman correlation with a constant vector is undefined (NaN) -- exclude
# degenerate genes explicitly upfront, with reporting, rather than let NaNs
# propagate silently into the thresholding step.
target_mat <- counts[targetscan_ids, , drop = FALSE]
gene_var <- apply(target_mat, 1, var)
n_zero_var <- sum(gene_var == 0 | is.na(gene_var))
cat(sprintf("Of %d candidate targets, %d have zero variance (excluded before correlation)\n",
            length(targetscan_ids), n_zero_var))
target_mat <- target_mat[gene_var > 0 & !is.na(gene_var), , drop = FALSE]
cat(sprintf("Proceeding with %d candidate targets for pairwise correlation\n", nrow(target_mat)))

# --- Validation check, now on real data -------------------------------------
set.seed(1)
val_idx <- sample(seq_len(nrow(target_mat)), min(10, nrow(target_mat)))
cpp_check <- spearman_full_cpp(target_mat[val_idx, , drop = FALSE])
r_check <- cor(t(target_mat[val_idx, , drop = FALSE]), method = "spearman")

max_diff <- max(abs(cpp_check - r_check))
cat(sprintf("\n*** VALIDATION: max |difference| between spearman_full_cpp() and R's cor() on %d sampled genes: %.4f ***\n",
            length(val_idx), max_diff))
if (max_diff > 0.01) {
  cat("*** WARNING: meaningful discrepancy found -- very likely the tie-handling\n",
      "difference flagged above. Consider fixing rank_vec() to average-rank ties\n",
      "before trusting the network built below, especially for lowly-expressed\n",
      "genes with many exact-zero values. ***\n")
} else {
  cat("Discrepancy negligible for this sample -- proceeding, but this was checked\n",
      "on only 10 genes; not a guarantee it holds for the full target set.\n")
}

# --- Full pairwise correlation matrix, via the C++ function -----------------
cat(sprintf("\nComputing all pairwise correlations among %d targets...\n", nrow(target_mat)))
cor_mat <- spearman_full_cpp(target_mat)
rownames(cor_mat) <- rownames(target_mat)
colnames(cor_mat) <- rownames(target_mat)

# --- Threshold: positive correlations > 0.4, excluding the diagonal --------
diag(cor_mat) <- 0  # remove self-correlation before thresholding
edge_mat <- cor_mat > CORR_THRESHOLD
n_edges <- sum(edge_mat) / 2  # matrix is symmetric, each edge counted twice
cat(sprintf("Edges with correlation > %.2f: %d\n", CORR_THRESHOLD, n_edges))

if (n_edges == 0) {
  stop(sprintf("No edges survive the > %.2f threshold. Check CORR_THRESHOLD, or whether the validation check above flagged a real discrepancy that needs fixing first.",
               CORR_THRESHOLD))
}

# --- Build graph, filter by minimum degree -----------------------------------
# Use only the upper triangle to avoid duplicate edges in the graph object.
edge_mat[lower.tri(edge_mat, diag = TRUE)] <- FALSE
edge_list <- which(edge_mat, arr.ind = TRUE)
edge_df <- data.frame(
  from = rownames(cor_mat)[edge_list[, 1]],
  to = rownames(cor_mat)[edge_list[, 2]]
)

g <- graph_from_data_frame(edge_df, directed = FALSE,
                            vertices = data.frame(name = rownames(cor_mat)))
node_degree <- degree(g)
cat(sprintf("Degree distribution: min=%d, median=%.0f, max=%d\n",
            min(node_degree), median(node_degree), max(node_degree)))

filtered_targets <- names(node_degree)[node_degree >= MIN_DEGREE]
cat(sprintf("Targets with degree >= %d: %d (of %d in the correlation network)\n",
            MIN_DEGREE, length(filtered_targets), nrow(target_mat)))

if (length(filtered_targets) < 5) {
  cat("*** WARNING: very few targets survived filtering -- results below may be unstable / not meaningful with so few genes. Consider lowering CORR_THRESHOLD or MIN_DEGREE. ***\n")
}

# --- Re-evaluate our method on the filtered target subset -------------------
compute_rho_for_genes <- function(gene_ids, label) {
  activity <- calculate_activity(counts = counts, targets = gene_ids)
  activity_df <- data.frame(activity = activity) |> rownames_to_column("Cell_ID")
  cor_df <- merge(activity_df, GFP_counts, by = "Cell_ID")
  cor_df <- cor_df[cor_df$eGFP < gfp_upper_cutoff, ]
  cor_df <- cor_df[is.finite(cor_df$activity) & is.finite(cor_df$eGFP), ]
  result <- cor.test(cor_df$activity, cor_df$eGFP, method = "spearman", exact = FALSE)
  cat(sprintf("%s: rho = %.3f, p = %.3e, n_targets = %d, n_cells = %d\n",
              label, result$estimate, result$p.value, length(gene_ids), nrow(cor_df)))
  list(cor_df = cor_df, cor = result)
}

cat("\n=== Results ===\n")
baseline_all_res <- compute_rho_for_genes(targetscan_ids, "Baseline (all TargetScan candidates, unfiltered)")
filtered_res <- compute_rho_for_genes(filtered_targets, sprintf("Co-expression-filtered (degree >= %d)", MIN_DEGREE))

N_TOP_TARGETS <- 200
targetscan_sorted <- targetscan_all[order(targetscan_all$Cumulative.weighted.context...score), ]
targets_top200 <- intersect(targetscan_sorted[1:N_TOP_TARGETS, ]$ensembl_gene_id, rownames(counts))
top200_res <- compute_rho_for_genes(targets_top200, "Baseline (TargetScan top 200, for reference)")

cat(sprintf("\n=== Summary ===\n"))
cat(sprintf("All TargetScan candidates (%d genes):        rho = %.3f\n",
            length(targetscan_ids), baseline_all_res$cor$estimate))
cat(sprintf("TargetScan top 200 (reference):                rho = %.3f\n", top200_res$cor$estimate))
cat(sprintf("Co-expression-filtered (%d genes, degree>=%d): rho = %.3f\n",
            length(filtered_targets), MIN_DEGREE, filtered_res$cor$estimate))

saveRDS(list(
  filtered_targets = filtered_targets,
  graph = g,
  node_degree = node_degree,
  cor_mat = cor_mat,
  baseline_all_res = baseline_all_res,
  filtered_res = filtered_res,
  top200_res = top200_res,
  corr_threshold = CORR_THRESHOLD,
  min_degree = MIN_DEGREE,
  cpp_validation_max_diff = max_diff
), file.path(OUT_DIR, "res_coexpression_filtered_targets.rds"))

cat(sprintf("\nSaved to %s\n", file.path(OUT_DIR, "res_coexpression_filtered_targets.rds")))
