library(tidyverse)
library(Matrix)
library(future)
library(furrr)
library(SingleCellExperiment)
library(scran)
library(scuttle)

# ---------------------------------------------------------------------------
# SCRAN, FIXED-N=100 DESIGN -- rebuilt for the new directory structure
# (/cfs/klemming/projects/snic/naiss2024-6-235/miRNA_footprint_new) and to
# match build_annotated_scran_dataset.R's exact scran recipe (cached
# protein-coding list, MAX_ZERO_PROP=0.75, dynamic
# GFP name resolution, log=FALSE linear-scale normcounts).
#
# CHANGES FROM THE PREVIOUS DRAFT, AND WHY:
#
# 1. Protein-coding list loaded from the CACHED biomart_protein_coding_ids.rds,
#    never a live biomaRt query. A live query would very likely FAIL on a
#    SLURM compute node (no outbound internet on compute nodes, typically
#    only login nodes have it) -- this isn't a style choice, it's a real
#    practical requirement for this to run via sbatch at all.
#
# 2. GFP ground truth now comes from cell_metadata.rds's GFP_normalized
#    column (scran-normalized, computed once on the full dataset), not a
#    standalone eGFP_counts.rds -- that file doesn't exist in the new
#    structure. GFP_normalized is the direct scran-world analog of the CPM
#    version's normalized ground truth.
#
# 3. GFP gene name resolved dynamically (same POSSIBLE_GFP_NAMES logic as
#    build_annotated_scran_dataset.R), not hardcoded to "eGFP" -- that
#    script explicitly flagged the hardcoded version as a latent risk.
#
# 4. Two exact-count stopifnot()s from the old draft (nrow(other_raw)==43960,
#    sum(ercc_idx)==39) are REMOVED, not carried over -- they were validated
#    against the OLD combined_counts_filtered_UMIs.rds. The regenerated
#    three-plate file may legitimately have different counts (built via
#    Reduce(intersect, ...) across plates). Replaced with informative
#    reporting instead of a hard crash on an expected difference.
#
# 5. FIXED n=100 cells only, per explicit decision -- no "all eligible
#    cells" arm, no auto-derived MATCHED_N. n=100 was chosen after n=50
#    repeatedly failed scran's clustering/pooling requirements -- 100
#    lands exactly on scran's own documented "we would want at least 100
#    cells per cluster" guidance, which may resolve this more cleanly than
#    continuing to force smaller pool sizes at n=50 did. Not yet confirmed
#    empirically -- treat as a stronger candidate, not a guaranteed fix.
# ---------------------------------------------------------------------------

BASE_DIR <- "~/Desktop/Projects/miRNA_footprint_new/"
PROCESSED_DIR <- file.path(BASE_DIR, "datasets/HEK_SS3/processed")
OUT_DIR <- file.path(BASE_DIR, "analysis/HEK_SS3/depth_subsampling_scran")
dir.create(file.path(OUT_DIR, "plots"), recursive = TRUE, showWarnings = FALSE)

# TODO: PLACEHOLDER -- calculate_activity() (called below) is not yet
# defined anywhere in this script. Utils.R's location in the new structure
# is unconfirmed (not visible in the directory tree shared so far). Update
# this path once confirmed -- the script will fail immediately at the first
# calculate_activity() call otherwise, not silently.
source(file.path(BASE_DIR, "scripts/Utils.R"))  # <-- UNVERIFIED PATH, FIX BEFORE RUNNING

QUICK_TEST_MODE <- FALSE   # TRUE = 2 depth levels x 3 iterations, to validate before the real sbatch job
N_WORKERS <- 8            # SET THIS to match your SLURM --cpus-per-task allocation

GFP_UPPER_PCT <- 0.99
N_TOP_TARGETS <- 200
MAX_ZERO_PROP <- 0.75          # matches build_annotated_scran_dataset.R
# SCRAN_CLUSTER_SEED removed -- was freezing scran's internal randomness
# identically across every iteration, hiding a real, large source of
# variance (see design note in compute_activity_gfp_cor_scran()). scran
# now inherits randomness from furrr's own per-task RNG stream instead.
POSSIBLE_GFP_NAMES <- c("eGFP", "GFP", "Gfp", "EGFP", "EGfp")

FIXED_N <- 100
N_ITER <- if (QUICK_TEST_MODE) 3 else 100
depth_grid <- if (QUICK_TEST_MODE) c(100, 100000) else c(1000,2000,3000,4000, 5000, 10000,20000, 50000, 100000)
cat("Depth levels to test:\n"); print(depth_grid)
cat(sprintf("N_ITER = %d, FIXED_N = %d cells, QUICK_TEST_MODE = %s\n", N_ITER, FIXED_N, QUICK_TEST_MODE))

use_scuttle <- requireNamespace("scuttle", quietly = TRUE)
cat(sprintf("Using scuttle::downsampleMatrix for exact thinning: %s\n", use_scuttle))

# --- Load raw (pre-normalization) counts -----------------------------------

raw <- readRDS(file.path(PROCESSED_DIR, "combined_counts_filtered_UMIs.rds"))
stopifnot(inherits(raw, "dgCMatrix"))

# --- Resolve GFP gene name dynamically (see design note 3) ------------------
egfp_gene_id <- NULL
for (name in POSSIBLE_GFP_NAMES) {
  if (name %in% rownames(raw)) { egfp_gene_id <- name; break }
}
if (is.null(egfp_gene_id)) {
  gfp_matches <- grep("GFP|gfp|Gfp", rownames(raw), value = TRUE)
  if (length(gfp_matches) == 0) stop("No GFP gene found in combined_counts_filtered_UMIs.rds!")
  egfp_gene_id <- gfp_matches[1]
}
cat(sprintf("Resolved GFP gene ID: %s\n", egfp_gene_id))

# Strip Ensembl version suffixes.
rn <- rownames(raw)
ensembl_idx <- grepl("^ENSG", rn)
rn[ensembl_idx] <- sub("\\.\\d+$", "", rn[ensembl_idx])
stopifnot(anyDuplicated(rn) == 0)
rownames(raw) <- rn

eGFP_idx <- which(rownames(raw) == egfp_gene_id)
stopifnot(length(eGFP_idx) == 1)
ercc_idx_full <- grepl("^ERCC", rownames(raw))
other_idx_full <- setdiff(seq_len(nrow(raw)), eGFP_idx)
other_idx_full <- other_idx_full[!grepl("^ERCC", rownames(raw)[other_idx_full])]

full_raw <- raw
other_raw <- raw[other_idx_full, , drop = FALSE]

cat(sprintf("Full matrix (denominator basis): %d genes x %d cells\n", nrow(full_raw), ncol(full_raw)))
cat(sprintf("Endogenous-only matrix: %d genes x %d cells (%d ERCC rows excluded)\n",
            nrow(other_raw), ncol(other_raw), sum(ercc_idx_full)))
# NOTE: no exact-count assertion here (see design note 4) -- report and move
# on, don't crash on a legitimately different count from the regenerated file.

# --- Protein-coding restriction, from the CACHED list (design note 1) ------
pc_cache_path <- file.path(PROCESSED_DIR, "biomart_protein_coding_ids.rds")
if (!file.exists(pc_cache_path)) {
  stop(sprintf("Expected cached protein-coding gene list at %s but it doesn't exist. ",
               "This script deliberately does NOT fall back to a live biomaRt query -- ",
               "that would very likely fail on a SLURM compute node with no internet access. ",
               "Run build_annotated_scran_dataset.R (or otherwise populate this cache) first.", pc_cache_path))
}
pc_query <- readRDS(pc_cache_path)
protein_coding_ids <- pc_query$ensembl_gene_id
other_raw_PC <- other_raw[rownames(other_raw) %in% protein_coding_ids, , drop = FALSE]
cat(sprintf("Endogenous genes restricted to protein-coding: %d (of %d)\n", nrow(other_raw_PC), nrow(other_raw)))

# --- GFP ground truth from cell_metadata.rds (design note 2) ---------------
cell_metadata <- readRDS(file.path(PROCESSED_DIR, "cell_metadata.rds"))
stopifnot(all(c("cell_id", "GFP_normalized") %in% colnames(cell_metadata)))
GFP_counts <- cell_metadata[!is.na(cell_metadata$GFP_normalized), c("cell_id", "GFP_normalized")]
colnames(GFP_counts) <- c("Cell_ID", "eGFP")
gfp_upper_cutoff <- quantile(GFP_counts$eGFP, GFP_UPPER_PCT)
cat(sprintf("Fixed GFP upper cutoff for QC exclusion only (scran-normalized, full dataset, %.0fth pct): %.3f\n",
            GFP_UPPER_PCT * 100, gfp_upper_cutoff))
cat(sprintf("GFP_normalized available for %d / %d cells (NA for cells excluded from the full-dataset scran run)\n",
            nrow(GFP_counts), nrow(cell_metadata)))

# --- Targets -----------------------------------------------------------------
targetscan_all <- readRDS(file.path(BASE_DIR, "resources/Targets__combined_124-3p_124-3p.2_506-3p.rds"))
targetscan_sorted <- targetscan_all[order(targetscan_all$Cumulative.weighted.context...score), ]
targets_top <- targetscan_sorted[1:N_TOP_TARGETS, ]$ensembl_gene_id

n_targets_present <- sum(targets_top %in% rownames(other_raw_PC))
cat(sprintf("Of top %d targets, %d present in protein-coding-restricted matrix (%d missing)\n",
            N_TOP_TARGETS, n_targets_present, N_TOP_TARGETS - n_targets_present))

# --- Original per-cell depth (eligibility + thinning probability) ---------
orig_depth <- Matrix::colSums(full_raw)
cat("Original per-cell depth (full library):\n"); print(summary(orig_depth))

n_eligible_per_depth <- sapply(depth_grid, function(d) sum(orig_depth >= d))
names(n_eligible_per_depth) <- depth_grid
cat("Eligible cells per depth level (pre-scan, before thinning):\n"); print(n_eligible_per_depth)
if (any(n_eligible_per_depth < FIXED_N)) {
  stop(sprintf("At least one depth level has fewer than FIXED_N=%d eligible cells (see table above). ",
               "Cannot proceed with a fixed-n=%d design at that depth -- either drop that depth level ",
               "or reduce FIXED_N.", FIXED_N, FIXED_N))
}

# --- Thinning function (unchanged mechanics from the CPM version) ---------
thin_matrix <- function(mat, target_depth) {
  cell_depths <- Matrix::colSums(mat)
  eligible <- which(cell_depths >= target_depth)
  if (length(eligible) == 0) return(NULL)
  mat_sub <- mat[, eligible, drop = FALSE]

  if (use_scuttle) {
    thinned <- scuttle::downsampleMatrix(mat_sub, prop = target_depth / cell_depths[eligible], bycol = TRUE)
  } else {
    thinned <- mat_sub
    probs_per_cell <- target_depth / cell_depths[eligible]
    col_of_entry <- rep(seq_len(ncol(mat_sub)), diff(mat_sub@p))
    thinned@x <- as.double(rbinom(n = length(mat_sub@x), size = mat_sub@x, prob = probs_per_cell[col_of_entry]))
    thinned <- Matrix::drop0(thinned)
  }
  list(matrix = thinned, eligible_cells = colnames(mat_sub))
}

# --- Correlation helper: scran, fixed-n only --------------------------------
compute_activity_gfp_cor_scran <- function(thinned_full, use_cells) {
  sub_mat <- thinned_full[, use_cells, drop = FALSE]
  new_totals <- Matrix::colSums(sub_mat)
  keep <- new_totals > 0
  sub_mat <- sub_mat[, keep, drop = FALSE]

  if (ncol(sub_mat) < 4) {
    return(list(rho = NA_real_, p_value = NA_real_, n_used = ncol(sub_mat), frac_zero_gfp = NA_real_,
                scran_failed = FALSE, n_genes_for_scran = 0, n_targets_used = NA_integer_))
  }
  cell_ids <- colnames(sub_mat)

  # Densified here deliberately: comparing a sparse dgCMatrix directly to 0
  # (pc_endog == 0) is rejected by newer Bioconductor SparseArray backends
  # (would silently force a huge dense result, so it's blocked instead --
  # confirmed by the actual error hit on first run, not assumed in advance).
  # Densifying is cheap at THIS scale (protein-coding genes x FIXED_N=100
  # cells, not the full ~44k x 753 matrix) -- unlike earlier in this
  # project, staying sparse doesn't matter here anymore.
  pc_endog <- as.matrix(sub_mat[rownames(sub_mat) %in% rownames(other_raw_PC), , drop = FALSE])
  zero_prop <- rowMeans(pc_endog == 0)
  pc_endog_filt <- pc_endog[zero_prop < MAX_ZERO_PROP, , drop = FALSE]

  if (nrow(pc_endog_filt) < 50) {   # gene-count floor: reasonable heuristic, not a verified scran requirement
    return(list(rho = NA_real_, p_value = NA_real_, n_used = 0, frac_zero_gfp = NA_real_,
                scran_failed = TRUE, n_genes_for_scran = nrow(pc_endog_filt),
                n_targets_used = sum(targets_top %in% rownames(pc_endog_filt))))
  }

  size_factors <- tryCatch({
    sce <- SingleCellExperiment(assays = list(counts = pc_endog_filt))
    # NO set.seed() here, deliberately -- removed after finding that
    # freezing this to a constant (previously SCRAN_CLUSTER_SEED=123 on
    # EVERY call) hid a real source of variance: computeSumFactors()'s
    # pooling/solving has a stochastic component, and results changed
    # substantially across different seeds. Freezing it meant our 20
    # per-depth iterations only varied in cell-sampling/thinning
    # randomness, silently excluding scran's OWN estimation noise from
    # what our boxplots report -- inconsistent with how every other noise
    # source in this project has been treated (always averaged over
    # explicitly, never hidden). future_map_dfr()'s furrr_options(seed=1312)
    # already gives each (depth, iteration) task its own well-separated,
    # reproducible RNG stream -- letting scran inherit from THAT, rather
    # than overriding it with one constant, means the 20-iteration spread
    # now honestly includes scran's own instability as part of the
    # reported uncertainty, which is arguably the more important thing
    # this analysis should be showing given how large that instability
    # turned out to be.
    # NO quickCluster() here, deliberately -- its default min.size=100
    # would consume the entire cell pool with nothing left over for
    # multiple clusters. Following scran author Aaron Lun's documented
    # guidance for small-n data: treat all cells as one pool via
    # computeSumFactors() directly.
    #
    # sizes= explicitly reduced from the default (c(20,40,60,80,100), max
    # 100) since a pool size equal to the WHOLE cell count leaves no room
    # for the sliding-window pooling the method relies on. c(20,30,40) was
    # chosen to keep max(sizes)=40 comfortably below n=100 with margin,
    # informed by (but not fully confirmed against) a "pool count should
    # be at least ~2x the largest pool size" heuristic seen in scran
    # discussions. HONESTY NOTE: the exact constraint that caused repeated
    # failures at n=50 was never fully pinned down against the currently
    # installed scran version -- if this STILL errors, run the constraint
    # down directly and interactively (see chat) rather than guess a
    # fourth value blind.
    sce <- computeSumFactors(sce, sizes = c(20, 30, 40))
    sf <- sizeFactors(sce)
    if (any(is.na(sf)) || any(sf <= 0)) stop("invalid (NA or non-positive) size factors")
    sf
  }, error = function(e) {
    cat(sprintf("  scran failed for this draw (%d cells, %d genes): %s\n",
                ncol(pc_endog_filt), nrow(pc_endog_filt), conditionMessage(e)))
    NULL
  })

  if (is.null(size_factors)) {
    return(list(rho = NA_real_, p_value = NA_real_, n_used = 0, frac_zero_gfp = NA_real_,
                scran_failed = TRUE, n_genes_for_scran = nrow(pc_endog_filt),
                n_targets_used = sum(targets_top %in% rownames(pc_endog_filt))))
  }
  names(size_factors) <- cell_ids

  scran_norm_endog <- t(t(pc_endog_filt) / size_factors)
  activity <- calculate_activity(counts = as.data.frame(scran_norm_endog), targets = targets_top)
  activity_df <- data.frame(activity = activity)
  n_targets_this_draw <- sum(targets_top %in% rownames(pc_endog_filt))

  eGFP_row_idx <- which(rownames(sub_mat) == egfp_gene_id)
  thinned_gfp_raw <- setNames(as.numeric(sub_mat[eGFP_row_idx, ]), cell_ids)
  frac_zero_gfp <- mean(thinned_gfp_raw == 0)
  thinned_gfp_scran <- thinned_gfp_raw / size_factors
  gfp_df <- data.frame(GFP = thinned_gfp_scran)

  cor_df <- merge(activity_df, gfp_df, by = 0)
  rownames(cor_df) <- cor_df$Row.names
  cor_df <- cor_df[, c("activity", "GFP")]

  orig_gfp_matched <- GFP_counts$eGFP[match(rownames(cor_df), GFP_counts$Cell_ID)]
  qc_pass <- !is.na(orig_gfp_matched) & orig_gfp_matched < gfp_upper_cutoff
  cor_df <- cor_df[qc_pass, ]
  cor_df <- cor_df[is.finite(cor_df$activity) & is.finite(cor_df$GFP), ]

  if (nrow(cor_df) < 4 || length(unique(cor_df$activity)) < 2 || length(unique(cor_df$GFP)) < 2) {
    return(list(rho = NA_real_, p_value = NA_real_, n_used = nrow(cor_df), frac_zero_gfp = frac_zero_gfp,
                scran_failed = FALSE, n_genes_for_scran = nrow(pc_endog_filt), n_targets_used = n_targets_this_draw))
  }

  cor <- cor.test(cor_df$activity, cor_df$GFP, method = "spearman", exact = FALSE)
  list(rho = unname(cor$estimate), p_value = cor$p.value, n_used = nrow(cor_df), frac_zero_gfp = frac_zero_gfp,
       scran_failed = FALSE, n_genes_for_scran = nrow(pc_endog_filt), n_targets_used = n_targets_this_draw)
}

# --- Depth-subsampling loop, fixed n=100 only -------------------------------

plan(multisession, workers = N_WORKERS)

depth_results <- future_map_dfr(
  depth_grid,
  function(target_depth) {
    map_dfr(1:N_ITER, function(iter) {
      th <- thin_matrix(full_raw, target_depth)
      if (is.null(th) || length(th$eligible_cells) < FIXED_N) {
        return(tibble(target_depth = target_depth, iteration = iter,
                       rho = NA_real_, p_value = NA_real_, n_used = 0, frac_zero_gfp = NA_real_,
                       scran_failed = NA, n_genes_for_scran = NA_integer_, n_targets_used = NA_integer_,
                       n_eligible_cells = if (is.null(th)) 0 else length(th$eligible_cells)))
      }

      thinned_full <- th$matrix
      fixed_cells <- sample(th$eligible_cells, size = FIXED_N, replace = FALSE)
      res <- compute_activity_gfp_cor_scran(thinned_full, fixed_cells)

      tibble(target_depth = target_depth, iteration = iter,
             rho = res$rho, p_value = res$p_value, n_used = res$n_used, frac_zero_gfp = res$frac_zero_gfp,
             scran_failed = res$scran_failed, n_genes_for_scran = res$n_genes_for_scran,
             n_targets_used = res$n_targets_used, n_eligible_cells = length(th$eligible_cells))
    })
  },
  .options = furrr_options(seed = 1312)
)

result_suffix <- if (QUICK_TEST_MODE) "_QUICKTEST" else ""
saveRDS(depth_results, file.path(OUT_DIR, sprintf("res_depth_subsampling_scran_n%d%s.rds", FIXED_N, result_suffix)))
depth_results <- readRDS("~/Desktop/Projects/miRNA_footprint_new/analysis/HEK_SS3/depth_subsampling_scran/res_depth_subsampling_scran_n100.rds")
# --- Failure / NA accounting (explicit) -------------------------------------
cat("\n=== scran failure rate (scran itself erroring, distinct from ordinary NA) ===\n")
depth_results |>
  group_by(target_depth) |>
  summarise(frac_scran_failed = mean(scran_failed, na.rm = TRUE),
            mean_n_genes_for_scran = mean(n_genes_for_scran, na.rm = TRUE),
            mean_n_targets_used = mean(n_targets_used, na.rm = TRUE),
            .groups = "drop") |>
  print()

n_na <- sum(is.na(depth_results$rho))
cat(sprintf("\nOverall: %d / %d draws were NA (%.2f%%)\n", n_na, nrow(depth_results), 100 * n_na / nrow(depth_results)))
if (n_na > 0) print(depth_results |> dplyr::filter(is.na(rho)) |> count(target_depth, name = "n_na"))

# --- Summarize + plot --------------------------------------------------------
summary_df <- depth_results |>
  group_by(target_depth) |>
  summarise(
    mean_rho = mean(rho, na.rm = TRUE),
    lower_rho = quantile(rho, 0.025, na.rm = TRUE),
    upper_rho = quantile(rho, 0.975, na.rm = TRUE),
    n_valid_draws = sum(!is.na(rho)),
    mean_n_targets_used = mean(n_targets_used, na.rm = TRUE),
    .groups = "drop"
  )
cat("\nSummary:\n"); print(summary_df)

y_range <- range(c(summary_df$lower_rho, summary_df$upper_rho, depth_results$rho), na.rm = TRUE)
y_limits <- c(min(0, floor(y_range[1] * 20) / 20), max(1, ceiling(y_range[2] * 20) / 20))
# Extend the top of the range slightly so the target-count labels have room
# above the highest box/whisker without getting clipped.
label_y <- y_limits[2] - 0.03 * diff(y_limits)
y_limits[2] <- y_limits[2] + 0.06 * diff(y_limits)

p_rho <- ggplot(depth_results, aes(x = factor(target_depth), y = rho)) +
  geom_boxplot(fill = "#8e7cc3", alpha = 0.7, outlier.shape = NA, width = 0.6) +
  geom_jitter(width = 0.15, alpha = 0.3, size = 0.8) +
  geom_text(data = summary_df, aes(x = factor(target_depth), y = label_y,
                                    label = sprintf("%.1f", mean_n_targets_used)),
            inherit.aes = FALSE, size = 3, vjust = 0) +
  theme_bw(base_size = 16) +
  theme(panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1)) +
  ylim(y_limits[1], y_limits[2]) +
  labs(x = "Target reads per cell (downsampled; evenly-spaced categories)", y = "Spearman rho",
       title = sprintf("Prediction accuracy vs. sequencing depth (scran, fixed n=%d cells)", FIXED_N))

print(p_rho)
ggsave(file.path(OUT_DIR, "plots", sprintf("depth_subsampling_rho_scran_n%d%s.pdf", FIXED_N, result_suffix)),
       p_rho, width = 8, height = 6, dpi = 600)

cat(sprintf("\nDone. Outputs written to %s\n", OUT_DIR))
if (QUICK_TEST_MODE) {
  cat("\n*** QUICK_TEST_MODE was TRUE -- this was a 2-depth, 3-iteration validation run only. ***\n",
      "*** Set QUICK_TEST_MODE <- FALSE for the real job before submitting via sbatch. ***\n")
}
