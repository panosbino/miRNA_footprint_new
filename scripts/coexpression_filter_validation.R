suppressMessages({
  library(dplyr); library(tibble); library(tidyr); library(ggplot2)
  library(Rcpp); library(RcppArmadillo)
})
# igraph is used through its namespace (igraph::...), not attached, so it
# cannot mask dplyr / base functions.
stopifnot(requireNamespace("igraph", quietly = TRUE))

# ===========================================================================
# Validation of the co-expression target filter (miR-124-3p, HEK293, scran)
#
# The filter: among TargetScan candidate targets, keep genes that have >= 3
# neighbours with Spearman correlation > 0.4 (computed across cells). Its
# first result (35 genes, rho 0.762 vs GFP) was computed with a C++ rank
# function that did not average tied ranks; that function is now fixed, and
# this script (1) proves the fix, then asks whether the result means anything:
#
#   [1] C++ correlation matches R's cor(method = "spearman"), including at
#       the level of which pairs cross the 0.4 threshold.
#   [2] Size-matched random null: does a co-expression-selected set of K
#       genes beat random sets of K genes from the SAME candidate pool, and
#       does it beat TargetScan's own top-K at equal size?
#   [3] Random-pool control: run the identical filter on expression-matched
#       NON-target genes. If those pools also give high-rho survivors, the
#       filter is picking up generic co-expression, not miR-124 targets.
#   [4] What distinguishes the selected genes from the other candidates
#       (TargetScan score, 8mer sites, 3'UTR length, expression level)?
#
# GFP is never used to select genes anywhere in this script; it is only used
# afterwards to score a gene set. Everything runs on all cells for network
# building, and on the usual GFP-QC cells (below the 99th percentile) for
# scoring, exactly as in the other comparisons.
# ===========================================================================

env_int <- function(name, default) as.integer(Sys.getenv(name, as.character(default)))

BASE_DIR      <- Sys.getenv("COEXP_BASE_DIR", "/cfs/klemming/projects/supr/naiss2024-6-235/miRNA_footprint_new")
PROCESSED_DIR <- file.path(BASE_DIR, "datasets/HEK_SS3/processed")
OUT_DIR       <- file.path(BASE_DIR, "analysis/HEK_SS3/comparisons")
TARGETS_PATH  <- file.path(BASE_DIR, "resources/human/Targets__combined_124-3p_124-3p.2_506-3p.rds")
UTILS_PATH    <- file.path(BASE_DIR, "scripts/HEK_SS3/Utils.R")          # defines calculate_activity()
CPP_PATH      <- file.path(BASE_DIR, "scripts/HEK_SS3/fullSpearman.cpp") # MUST be the corrected version
UTR_PATH      <- file.path(BASE_DIR, "tools/miReact/seqs/hs.utr3.seqs.rds")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

GFP_UPPER_PCT  <- 0.99
CORR_THRESHOLD <- 0.4
MIN_DEGREE     <- 3
N_TOP_TARGETS  <- 200                       # the standard TargetScan top-200 used in every other comparison
N_NULL         <- env_int("COEXP_N_NULL", 1000)   # random K-subsets of the candidate pool   [item 2]
N_POOLS        <- env_int("COEXP_N_POOLS", 100)   # expression-matched random pools           [item 3]
N_EXPR_BINS    <- 30                        # expression-matching resolution for random pools
MIN_SURV       <- 5                         # fewer survivors than this => rho not computed for that pool
COR_TOL        <- 1e-6
SEED           <- 2026
set.seed(SEED)

source(UTILS_PATH)

# The C++ function is only a speed-up: R's own cor(method = "spearman") gives the
# same numbers (max difference ~2e-16) and costs ~2.5x more time here (about 3 min
# instead of 1 for all random pools). So if the C++ cannot be compiled (e.g. on a
# Mac without the gfortran 14.2 runtime that CRAN's R build links against), fall
# back to R rather than stop. Force the fallback with COEXP_FORCE_R_CORRELATION=1.
USE_CPP <- FALSE
if (Sys.getenv("COEXP_FORCE_R_CORRELATION", "0") != "1") {
  USE_CPP <- tryCatch({ sourceCpp(CPP_PATH); TRUE }, error = function(e) {
    cat("\n*** Could not compile/link ", CPP_PATH, " (", conditionMessage(e), ").\n",
        "    Falling back to R's cor(method = 'spearman'): identical results, a bit slower.\n",
        "    To use the C++ on macOS, install the gfortran 14.2 package from https://mac.r-project.org/tools/ ***\n\n", sep = "")
    FALSE
  })
}
if (!USE_CPP) spearman_full_cpp <- function(X) cor(t(X), method = "spearman")   # genes x cells in, genes x genes out

# --- Data ----------------------------------------------------------------------
counts <- as.matrix(readRDS(file.path(PROCESSED_DIR, "scran_normalized_linear.rds")))
cat(sprintf("Counts: %d genes x %d cells (scran-normalized)\n", nrow(counts), ncol(counts)))

cell_metadata <- readRDS(file.path(PROCESSED_DIR, "cell_metadata.rds"))
GFP_counts <- cell_metadata[!is.na(cell_metadata$GFP_normalized), c("cell_id", "GFP_normalized")]
colnames(GFP_counts) <- c("Cell_ID", "eGFP")
gfp_upper_cutoff <- quantile(GFP_counts$eGFP, GFP_UPPER_PCT)
qc <- GFP_counts[GFP_counts$eGFP < gfp_upper_cutoff & GFP_counts$Cell_ID %in% colnames(counts), ]
cells_qc <- qc$Cell_ID
gfp_qc   <- qc$eGFP
cat(sprintf("Scoring cells after GFP QC: %d\n", length(cells_qc)))

targetscan_all <- readRDS(TARGETS_PATH)
stopifnot(all(c("ensembl_gene_id", "Cumulative.weighted.context...score") %in% colnames(targetscan_all)))

cand_ids <- unique(intersect(targetscan_all$ensembl_gene_id, rownames(counts)))
cand_var <- apply(counts[cand_ids, , drop = FALSE], 1, var)
cand_ids <- cand_ids[is.finite(cand_var) & cand_var > 0]
cat(sprintf("TargetScan candidates: %d in file, %d present in the scran matrix (non-zero variance)\n",
            length(unique(targetscan_all$ensembl_gene_id)), length(cand_ids)))
cat(sprintf("Fraction of exact zeros in the candidate matrix: %.3f (ties are common, which is why rank handling matters)\n",
            mean(counts[cand_ids, ] == 0)))

# --- Scoring functions -----------------------------------------------------------
# Reference route: exactly the established pipeline (calculate_activity ->
# merge with GFP -> QC -> Spearman).
rho_ref <- function(genes) {
  act <- calculate_activity(counts = counts, targets = genes)
  adf <- data.frame(Cell_ID = names(act), activity = as.numeric(act))
  m <- merge(adf, GFP_counts, by = "Cell_ID")
  m <- m[m$eGFP < gfp_upper_cutoff & is.finite(m$activity) & is.finite(m$eGFP), ]
  unname(cor.test(m$activity, m$eGFP, method = "spearman", exact = FALSE)$estimate)
}
# Fast route, for the thousands of null draws. activity = -log2(S / mean(S)) + 1
# is strictly DEcreasing in the target sum S, so Spearman(activity, GFP) is
# exactly -Spearman(S, GFP): the normalisation cannot change a rank
# correlation. This is verified against rho_ref below before it is used.
rho_fast <- function(genes) {
  genes <- intersect(genes, rownames(counts))
  if (length(genes) < 2) return(NA_real_)
  s <- colSums(counts[genes, cells_qc, drop = FALSE])
  ok <- is.finite(s) & s > 0            # mirrors dropping non-finite activity
  if (sum(ok) < 10) return(NA_real_)
  -suppressWarnings(cor(s[ok], gfp_qc[ok], method = "spearman"))
}

# --- The filter, as one reusable function ------------------------------------------
coexpr_select <- function(mat, thr = CORR_THRESHOLD, min_deg = MIN_DEGREE, keep_matrices = FALSE) {
  C <- spearman_full_cpp(mat)
  diag(C) <- 0
  A <- C > thr
  deg  <- rowSums(A)                    # number of neighbours above the threshold
  wdeg <- rowSums(C * A)                # summed correlation to those neighbours (tie-breaker)
  names(deg) <- names(wdeg) <- rownames(mat)
  out <- list(keep = names(deg)[deg >= min_deg], degree = deg, wdegree = wdeg, n_edges = sum(A) / 2)
  if (keep_matrices) { rownames(A) <- colnames(A) <- rownames(mat); out$adj <- A; out$cor <- C }
  out
}
top_k_by_degree <- function(sel, k) names(sel$degree)[order(-sel$degree, -sel$wdegree)][seq_len(k)]

# ===========================================================================
# [1] C++ VALIDATION
# ===========================================================================
cat("\n================ [1] C++ correlation validation ================\n")
val_ids <- sample(cand_ids, min(300, length(cand_ids)))   # drawn in both modes so the random stream is identical
if (USE_CPP) {
# (a) known-answer case with ties
Y <- rbind(c(1, 1, 1, 2, 3, 3), c(6, 5, 4, 3, 2, 1))
ka_cpp <- spearman_full_cpp(Y)[1, 2]; ka_r <- cor(Y[1, ], Y[2, ], method = "spearman")
cat(sprintf("known-answer case with ties: C++ %.10f vs R %.10f\n", ka_cpp, ka_r))
# (b) real candidate genes, including which pairs cross the edge threshold
Xv <- counts[val_ids, , drop = FALSE]
Cc <- spearman_full_cpp(Xv); Cr <- cor(t(Xv), method = "spearman")
ut <- upper.tri(Cc)
max_diff <- max(abs(Cc - Cr))
n_discordant <- sum((Cc[ut] > CORR_THRESHOLD) != (Cr[ut] > CORR_THRESHOLD))
cat(sprintf("%d real candidate genes: max |C++ - R| = %.2e; pairs on opposite sides of the %.1f threshold: %d of %d\n",
            length(val_ids), max_diff, CORR_THRESHOLD, n_discordant, sum(ut)))
if (abs(ka_cpp - ka_r) > COR_TOL || max_diff > COR_TOL || n_discordant > 0) {
  stop("spearman_full_cpp() does not reproduce R's cor(method = 'spearman'). Make sure ", CPP_PATH,
       " is the CORRECTED fullSpearman.cpp (average-rank ties); the original version fails exactly this check.")
}
cat("PASSED: the C++ function now agrees with R's Spearman, including at the edge threshold.\n")
} else {
  max_diff <- NA_real_; n_discordant <- NA_integer_
  cat("C++ not in use: all correlations below are computed with R's cor(method = 'spearman'), the reference\n",
      "implementation, so there is nothing to validate in this step.\n", sep = "")
}

# ===========================================================================
# OBSERVED RESULT (with the corrected function)
# ===========================================================================
cat("\n================ Observed co-expression filter ================\n")
cand_mat <- counts[cand_ids, , drop = FALSE]
obs <- coexpr_select(cand_mat, keep_matrices = TRUE)
K <- length(obs$keep)
cat(sprintf("Candidates: %d | edges (rho > %.1f): %d | genes with degree >= %d: %d\n",
            length(cand_ids), CORR_THRESHOLD, obs$n_edges, MIN_DEGREE, K))
cat(sprintf("Degree distribution: min=%d, median=%.0f, max=%d\n", min(obs$degree), median(obs$degree), max(obs$degree)))
if (K < MIN_SURV) stop("Fewer than ", MIN_SURV, " genes survive the filter with the corrected correlations; nothing further to validate.")

# the graph, as requested originally; its degrees must equal the matrix-based ones
g <- igraph::graph_from_adjacency_matrix(1 * obs$adj, mode = "undirected", diag = FALSE)
stopifnot(all(igraph::degree(g)[names(obs$degree)] == obs$degree))

# score-route cross-check: fast route must reproduce the established pipeline
ts_sorted <- targetscan_all[order(targetscan_all$Cumulative.weighted.context...score), ]
top200_genes <- intersect(ts_sorted[seq_len(N_TOP_TARGETS), ]$ensembl_gene_id, rownames(counts))
check_sets <- list(selected = obs$keep, all_candidates = cand_ids, top200 = top200_genes,
                   random_a = sample(cand_ids, K), random_b = sample(cand_ids, K))
chk <- sapply(check_sets, function(s) c(ref = rho_ref(s), fast = rho_fast(s)))
cat("\nScoring cross-check (established pipeline vs fast route used for the null draws):\n"); print(round(t(chk), 6))
if (max(abs(chk["ref", ] - chk["fast", ])) > 1e-8) {
  stop("Fast scoring route does not reproduce calculate_activity()-based rho. Do not use the null results; ",
       "check that Utils.R's calculate_activity() matches -log2(colSums/mean)+1.")
}
cat("PASSED: both routes agree.\n")

rho_obs        <- rho_ref(obs$keep)
rho_all        <- rho_ref(cand_ids)
rho_top200     <- rho_ref(top200_genes)
ts_score_gene  <- tapply(targetscan_all$Cumulative.weighted.context...score, targetscan_all$ensembl_gene_id, min, na.rm = TRUE)
ts_topK        <- names(sort(ts_score_gene[cand_ids]))[seq_len(K)]      # TargetScan's own top-K at the SAME size
rho_ts_topK    <- rho_ref(ts_topK)
cat(sprintf("\nrho vs GFP:\n  co-expression-selected (K = %d genes)      %.3f\n  TargetScan top-%d, same size                 %.3f\n  TargetScan top 200 (%d genes present)        %.3f\n  all %d candidates                          %.3f\n",
            K, rho_obs, K, rho_ts_topK, length(top200_genes), rho_top200, length(cand_ids), rho_all))

# ===========================================================================
# [2] SIZE-MATCHED RANDOM NULL
# ===========================================================================
cat(sprintf("\n================ [2] Random subsets of the candidate pool (K = %d, n = %d) ================\n", K, N_NULL))
null_rho <- replicate(N_NULL, rho_fast(sample(cand_ids, K)))
p_null <- (1 + sum(null_rho >= rho_obs, na.rm = TRUE)) / (1 + sum(!is.na(null_rho)))
cat(sprintf("random K-subsets: median %.3f | 95th pct %.3f | 99th pct %.3f | max %.3f\n",
            median(null_rho, na.rm = TRUE), quantile(null_rho, 0.95, na.rm = TRUE),
            quantile(null_rho, 0.99, na.rm = TRUE), max(null_rho, na.rm = TRUE)))
cat(sprintf("observed %.3f exceeds %.1f%% of random subsets; empirical p = %.4f (one-sided)\n",
            rho_obs, 100 * mean(null_rho < rho_obs, na.rm = TRUE), p_null))
cat(sprintf("Equal-size comparison with TargetScan's own top-%d: co-expression %.3f vs TargetScan %.3f (difference %+.3f); the two sets share %d of %d genes\n",
            K, rho_obs, rho_ts_topK, rho_obs - rho_ts_topK, length(intersect(obs$keep, ts_topK)), K))
cat("(if the two sets share most genes, the filter is largely re-finding TargetScan's strongest candidates and the difference says little)\n")
cat(sprintf("Smallest empirical p attainable with %d draws: %.4f\n", N_NULL, 1 / (1 + N_NULL)))

# ===========================================================================
# [3] RANDOM-POOL CONTROL (expression-matched non-target genes)
# ===========================================================================
cat(sprintf("\n================ [3] Random-pool control (%d expression-matched pools) ================\n", N_POOLS))
all_var <- apply(counts, 1, var)
mean_expr_all <- rowMeans(counts)
det_all <- rowMeans(counts > 0)
bin_all <- cut(rank(mean_expr_all, ties.method = "average") / length(mean_expr_all),
               breaks = seq(0, 1, length.out = N_EXPR_BINS + 1), include.lowest = TRUE, labels = FALSE)
names(bin_all) <- rownames(counts)

nontarget_ids <- setdiff(rownames(counts), targetscan_all$ensembl_gene_id)     # excludes ALL TargetScan candidates
nontarget_ids <- nontarget_ids[is.finite(all_var[nontarget_ids]) & all_var[nontarget_ids] > 0]
need <- table(factor(bin_all[cand_ids], levels = seq_len(N_EXPR_BINS)))
pool_by_bin <- split(nontarget_ids, factor(bin_all[nontarget_ids], levels = seq_len(N_EXPR_BINS)))
shortfall <- sum(pmax(as.integer(need) - lengths(pool_by_bin), 0))
cat(sprintf("Non-target gene pool: %d genes. Each random pool mirrors the candidates' expression distribution (%d bins).\n",
            length(nontarget_ids), N_EXPR_BINS))
if (shortfall > 0) cat(sprintf("NOTE: %d candidate slots have no available non-target gene in their bin; those pools will be slightly smaller.\n", shortfall))

draw_matched_pool <- function() {
  unlist(lapply(seq_len(N_EXPR_BINS), function(b) {
    n_need <- as.integer(need[b]); avail <- pool_by_bin[[b]]
    if (n_need == 0 || length(avail) == 0) return(NULL)
    sample(avail, min(n_need, length(avail)))
  }), use.names = FALSE)
}

obs_topK_set <- top_k_by_degree(obs, K)     # identical to obs$keep, since survivors outrank all other genes
rho_obs_topK <- rho_fast(obs_topK_set)

pool_res <- vector("list", N_POOLS)
for (b in seq_len(N_POOLS)) {
  pool <- draw_matched_pool()
  sel <- coexpr_select(counts[pool, , drop = FALSE])
  # All pool-level quantities are computed as plain variables BEFORE building the
  # row. (An earlier version named a tibble column `pool`, which then shadowed this
  # gene vector inside the later arguments of the same tibble() call.)
  pool_n_genes <- length(pool)
  pool_med_expr <- median(mean_expr_all[pool])
  pool_med_det <- median(det_all[pool])
  pool_rho_surv <- if (length(sel$keep) >= MIN_SURV) rho_fast(sel$keep) else NA_real_
  pool_rho_topK <- rho_fast(top_k_by_degree(sel, K))
  pool_res[[b]] <- tibble(pool_id = b, n_genes = pool_n_genes, n_edges = sel$n_edges,
                          n_survivors = length(sel$keep), rho_survivors = pool_rho_surv,
                          rho_topK = pool_rho_topK, med_expr = pool_med_expr, med_det = pool_med_det)
  if (b %% 10 == 0) cat(sprintf("  ... %d/%d pools done\n", b, N_POOLS))
}
pool_res <- bind_rows(pool_res)

cat(sprintf("\nBalance check (candidates vs matched pools): median mean-expression %.3f vs %.3f | median fraction detected %.3f vs %.3f | pool size %d vs %d candidates\n",
            median(mean_expr_all[cand_ids]), mean(pool_res$med_expr), median(det_all[cand_ids]), mean(pool_res$med_det),
            round(mean(pool_res$n_genes)), length(cand_ids)))
bal_rel <- abs(mean(pool_res$med_expr) - median(mean_expr_all[cand_ids])) / median(mean_expr_all[cand_ids])
if (bal_rel > 0.15) cat(sprintf("*** WARNING: matched pools differ from the candidates' median expression by %.0f%% -- the control is not well matched; do not interpret it. ***\n", 100 * bal_rel))
cat(sprintf("Edges (rho > %.1f):   observed candidates %d | random pools median %.0f (range %d-%d)\n",
            CORR_THRESHOLD, obs$n_edges, median(pool_res$n_edges), min(pool_res$n_edges), max(pool_res$n_edges)))
cat(sprintf("Survivors (degree >= %d): observed %d | random pools median %.0f (range %d-%d); %d of %d pools have >= %d\n",
            MIN_DEGREE, K, median(pool_res$n_survivors), min(pool_res$n_survivors), max(pool_res$n_survivors),
            sum(pool_res$n_survivors >= MIN_SURV), N_POOLS, MIN_SURV))
cat(sprintf("rho of the top-%d genes by degree (size-matched): observed %.3f | random pools median %.3f, range %.3f to %.3f\n",
            K, rho_obs_topK, median(pool_res$rho_topK, na.rm = TRUE),
            min(pool_res$rho_topK, na.rm = TRUE), max(pool_res$rho_topK, na.rm = TRUE)))
cat("(the sign matters: strongly co-expressed random genes can correlate with the activity score in EITHER direction,\n",
    " depending on how the shared factor relates to GFP; a large value of either sign means the filter finds generic structure)\n", sep = "")
if (any(!is.na(pool_res$rho_survivors))) {
  cat(sprintf("rho of the actual survivors in pools that have >= %d: median %.3f (n = %d pools)\n",
              MIN_SURV, median(pool_res$rho_survivors, na.rm = TRUE), sum(!is.na(pool_res$rho_survivors))))
}
p_pool_rho   <- (1 + sum(pool_res$rho_topK >= rho_obs_topK, na.rm = TRUE)) / (1 + sum(!is.na(pool_res$rho_topK)))
p_pool_edges <- (1 + sum(pool_res$n_edges >= obs$n_edges)) / (1 + N_POOLS)
cat(sprintf("empirical p: rho of top-%d by degree >= observed in random pools: %.4f | edges >= observed: %.4f   (smallest attainable with %d pools: %.4f)\n",
            K, p_pool_rho, p_pool_edges, N_POOLS, 1 / (1 + N_POOLS)))

# ===========================================================================
# [4] WHAT DISTINGUISHES THE SELECTED GENES FROM THE OTHER CANDIDATES?
# ===========================================================================
cat("\n================ [4] Selected genes vs other candidates ================\n")
feat <- tibble(ensembl_gene_id = cand_ids,
               mean_expr = rowMeans(cand_mat), detection = rowMeans(cand_mat > 0),
               selected = cand_ids %in% obs$keep)
feat$context_score <- as.numeric(ts_score_gene[cand_ids])

mer8_cols <- grep("8mer", colnames(targetscan_all), value = TRUE)
if (length(mer8_cols) > 0) {
  tmp <- targetscan_all[targetscan_all$ensembl_gene_id %in% cand_ids, ]
  M <- do.call(cbind, lapply(mer8_cols, function(cn) suppressWarnings(as.numeric(tmp[[cn]]))))
  tmp$mer8 <- rowSums(M, na.rm = TRUE)
  mer8_gene <- tapply(tmp$mer8, tmp$ensembl_gene_id, max)
  feat$mer8_sites <- as.numeric(mer8_gene[cand_ids])
  cat(sprintf("8mer site count = sum of: %s (max over rows per gene)\n", paste(mer8_cols, collapse = ", ")))
} else cat("No column containing '8mer' found in the TargetScan table; 8mer comparison skipped.\n")

if (file.exists(UTR_PATH)) {
  u <- readRDS(UTR_PATH)
  if (all(c("gid", "nchar") %in% colnames(u))) {
    utr_gene <- tapply(u$nchar, u$gid, max)
    feat$utr_length <- as.numeric(utr_gene[cand_ids])
    cat(sprintf("3'UTR length from %s (longest transcript per gene): available for %d of %d candidates\n",
                basename(UTR_PATH), sum(!is.na(feat$utr_length)), nrow(feat)))
  } else cat("UTR file lacks gid/nchar columns; UTR comparison skipped.\n")
} else cat("UTR file not found; UTR comparison skipped.\n")

compare_feature <- function(col, label) {
  if (!col %in% colnames(feat)) return(NULL)
  x <- feat[[col]][feat$selected]; y <- feat[[col]][!feat$selected]
  if (sum(is.finite(x)) < 3 || sum(is.finite(y)) < 3) return(NULL)
  tibble(feature = label, median_selected = median(x, na.rm = TRUE), median_other = median(y, na.rm = TRUE),
         wilcoxon_p = suppressWarnings(wilcox.test(x, y)$p.value))
}
feat_cmp <- bind_rows(
  compare_feature("context_score", "TargetScan context score (more negative = stronger)"),
  compare_feature("mer8_sites", "8mer sites"),
  compare_feature("utr_length", "3'UTR length (nt)"),
  compare_feature("mean_expr", "mean expression (scran)"),
  compare_feature("detection", "fraction of cells detected"))
cat("\nMedian per group (selected vs the other candidates), two-sided Wilcoxon:\n"); print(as.data.frame(feat_cmp), digits = 4, row.names = FALSE)

feat$ts_rank <- rank(feat$context_score, ties.method = "first")           # 1 = strongest predicted
enrich <- bind_rows(lapply(c(100, 200, 500), function(X) {
  top <- feat$ts_rank <= X; k <- sum(feat$selected); m <- sum(top); x <- sum(top & feat$selected)
  tibble(top_X_by_context_score = X, selected_in_top = x, expected_by_chance = k * m / nrow(feat),
         hypergeometric_p = phyper(x - 1, m, nrow(feat) - m, k, lower.tail = FALSE))
}))
cat(sprintf("\nOverlap of the %d selected genes with TargetScan's strongest candidates (ranked among the %d present):\n", K, nrow(feat)))
print(as.data.frame(enrich), digits = 4, row.names = FALSE)

# ===========================================================================
# PLOTS
# ===========================================================================
refs <- tibble(label = c(sprintf("co-expression selected (K=%d)", K), sprintf("TargetScan top-%d (same size)", K),
                         sprintf("TargetScan top 200 (%d genes)", length(top200_genes)),
                         sprintf("all %d candidates", length(cand_ids))),
               value = c(rho_obs, rho_ts_topK, rho_top200, rho_all))
p1 <- ggplot(tibble(rho = null_rho[!is.na(null_rho)]), aes(rho)) +
  geom_histogram(bins = 40, fill = "grey75", color = "grey40") +
  geom_vline(data = refs, aes(xintercept = value, color = label, linetype = label), linewidth = 1.1) +
  theme_bw(base_size = 14) + theme(legend.position = "bottom", legend.direction = "vertical", legend.title = element_blank()) +
  labs(x = "Spearman rho vs GFP", y = "random subsets",
       title = sprintf("Random %d-gene subsets of the candidate pool (n = %d)", K, N_NULL),
       subtitle = sprintf("empirical p for the co-expression-selected set = %.4f", p_null))
ggsave(file.path(OUT_DIR, "coexp_validation_random_subsets.pdf"), p1, width = 8, height = 7)

pool_long <- pool_res %>%
  transmute(`rho of top-K genes by degree` = rho_topK, `edges above threshold` = n_edges, `survivors (degree >= min)` = n_survivors) %>%
  pivot_longer(everything(), names_to = "metric", values_to = "value")
obs_lines <- tibble(metric = c("rho of top-K genes by degree", "edges above threshold", "survivors (degree >= min)"),
                    value = c(rho_obs_topK, obs$n_edges, K))
p2 <- ggplot(pool_long, aes(value)) +
  geom_histogram(bins = 25, fill = "grey75", color = "grey40", na.rm = TRUE) +
  geom_vline(data = obs_lines, aes(xintercept = value), color = "red", linewidth = 1.1) +
  facet_wrap(~metric, scales = "free", ncol = 1) + theme_bw(base_size = 14) +
  labs(x = NULL, y = "random pools", title = "Expression-matched random pools",
       subtitle = sprintf("n = %d pools; red line = observed candidate pool", N_POOLS))
ggsave(file.path(OUT_DIR, "coexp_validation_random_pools.pdf"), p2, width = 7, height = 9)

feat_long <- feat %>%
  transmute(selected = ifelse(selected, "selected", "other candidates"),
            `context score` = context_score,
            `log10 mean expression` = log10(mean_expr + 1e-6),
            `fraction detected` = detection,
            `log10 3'UTR length` = if ("utr_length" %in% colnames(feat)) log10(utr_length) else NA_real_,
            `8mer sites` = if ("mer8_sites" %in% colnames(feat)) mer8_sites else NA_real_) %>%
  pivot_longer(-selected, names_to = "feature", values_to = "value") %>% filter(is.finite(value))
p3 <- ggplot(feat_long, aes(selected, value, fill = selected)) +
  geom_boxplot(outlier.alpha = 0.3, show.legend = FALSE) + facet_wrap(~feature, scales = "free_y") +
  theme_bw(base_size = 13) + labs(x = NULL, y = NULL, title = "Selected genes vs other candidates")
ggsave(file.path(OUT_DIR, "coexp_validation_features.pdf"), p3, width = 9, height = 6)

saveRDS(list(params = list(corr_threshold = CORR_THRESHOLD, min_degree = MIN_DEGREE, K = K, n_null = N_NULL, n_pools = N_POOLS),
             selected = obs$keep, degree = obs$degree, n_edges = obs$n_edges, graph = g,
             rho = c(selected = rho_obs, ts_topK = rho_ts_topK, top200 = rho_top200, all_candidates = rho_all),
             null_rho = null_rho, p_null = p_null, pool_res = pool_res,
             p_pool_rho = p_pool_rho, p_pool_edges = p_pool_edges,
             features = feat, feature_comparison = feat_cmp, enrichment = enrich,
             cpp_validation = list(max_diff = max_diff, discordant_pairs = n_discordant)),
        file.path(OUT_DIR, "res_coexpression_validation.rds"))

cat("\n================ How to read this ================\n")
cat("[2] If the observed rho sits inside the bulk of the random-subset histogram, selecting by co-expression adds\n",
    "    nothing beyond picking candidates at random. If it is above nearly all draws, selection does something --\n",
    "    but not necessarily what we hope. Also compare with TargetScan's own top-K at the same size.\n",
    "[3] If expression-matched random pools often give top-K genes with rho comparable to the observed one, the filter\n",
    "    is mostly finding generic co-expression (e.g. a shared cell-state axis), not miR-124 targets. If random pools\n",
    "    have far fewer edges and clearly lower rho, the structure is specific to the candidate list. LIMIT OF THIS\n",
    "    CONTROL: it cannot tell miRNA-driven co-expression apart from co-expression due to other shared properties of\n",
    "    TargetScan candidates (similar 3'UTR length, tissue-specific modules); that needs another miRNA's candidate list.\n",
    "[4] Selected genes that look like ordinary strong candidates (strong score, more 8mer sites) support a target\n",
    "    interpretation; selected genes that are mainly highly expressed or long-UTR genes point to a technical or\n",
    "    compositional explanation instead.\n", sep = "")
cat(sprintf("\nSaved plots and res_coexpression_validation.rds to %s\n", OUT_DIR))
