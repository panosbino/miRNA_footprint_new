suppressMessages({
  library(dplyr); library(tibble); library(tidyr); library(ggplot2)
  library(Rcpp); library(RcppArmadillo)
})
# igraph is used through its namespace (igraph::...), not attached, so it cannot mask dplyr / base functions.
stopifnot(requireNamespace("igraph", quietly = TRUE))

# ===========================================================================
# Validation of the co-expression target filter, scored against GFP
# FLUORESCENCE (FACS index sort) on Dox-INDUCED cells only.
#
# What changed relative to coexpression_filter_validation.R:
#   * ground truth: GFP fluorescence (FITC-A) instead of GFP mRNA
#   * scoring cells: Dox-induced cells only (Dox_Group == "Dox Induced")
#   * the co-expression NETWORK can be built on two different cell sets, and
#     the script runs both by default so they can be compared:
#       "induced" : network from induced cells only. The induced-vs-uninduced
#                   step (and Dox exposure itself) cannot drive the selection.
#       "all"     : network from all cells = the ORIGINAL filter. Its selection
#                   is dominated by the induction step, and Dox changes gene
#                   expression independently of miR-124, so it may pick up
#                   Dox-responsive genes with no within-induced signal.
#   Both are scored identically, so differences are due to selection alone.
#
# The filter: among TargetScan candidates, keep genes with >= MIN_DEGREE
# neighbours at Spearman rho > CORR_THRESHOLD (computed across the network
# cells). Fluorescence is never used to select genes; only to score a set.
#
# Checks, per network variant:
#   [2] random K-subsets of the candidate pool (and TargetScan's own top-K)
#   [3] the identical filter on expression-matched NON-target pools
#   [4] what distinguishes the selected genes from the other candidates
# plus, once: [1] C++ vs R Spearman agreement, a scoring-route cross-check, and
# validation of the fluorescence-to-cell join (same code as
# fluorescence_accuracy_comparison.R).
# ===========================================================================

env_int <- function(name, default) as.integer(Sys.getenv(name, as.character(default)))
env_num <- function(name, default) as.numeric(Sys.getenv(name, as.character(default)))

BASE_DIR      <- Sys.getenv("COEXP_BASE_DIR", Sys.getenv("HEK_BASE_DIR", "~/Desktop/Projects/miRNA_footprint_new"))
RAW_DIR       <- file.path(BASE_DIR, "datasets/HEK_SS3/raw")
PROCESSED_DIR <- file.path(BASE_DIR, "datasets/HEK_SS3/processed")
OUT_DIR       <- file.path(BASE_DIR, "analysis/HEK_SS3/comparisons")
TARGETS_PATH  <- file.path(BASE_DIR, "resources/human/Targets__combined_124-3p_124-3p.2_506-3p.rds")
UTILS_PATH    <- file.path(BASE_DIR, "scripts/Utils.R")          # defines calculate_activity()
CPP_PATH      <- file.path(BASE_DIR, "scripts/fullSpearman.cpp") # MUST be the corrected version
UTR_PATH      <- file.path(BASE_DIR, "tools/miReact/seqs/hs.utr3.seqs.rds")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

SEQ_PREFIX      <- "P32156"                 # prefix of the per-plate sequencing files
PLATES          <- c("101", "103", "105")
INDUCED_GROUPS  <- "Dox Induced"            # 0.01, 0.1 and 1 Dox; excludes "Control", "0 Dox", "Unknown"
MIN_MATCH_FRAC  <- 0.5
MIN_RHO_SANITY  <- 0.2

CORR_THRESHOLD <- env_num("COEXP_CORR_THRESHOLD", 0.4)
MIN_DEGREE     <- env_int("COEXP_MIN_DEGREE", 3)
N_TOP_TARGETS  <- 200                               # the standard TargetScan top-200 used in every other comparison
N_NULL         <- env_int("COEXP_N_NULL", 1000)     # random K-subsets of the candidate pool   [2]
N_POOLS        <- env_int("COEXP_N_POOLS", 100)     # expression-matched random pools          [3]
N_EXPR_BINS    <- 30
MIN_SURV       <- 5                                 # fewer survivors than this => nothing to validate for that variant
COR_TOL        <- 1e-6
SEED           <- 2026
VARIANTS       <- trimws(strsplit(Sys.getenv("COEXP_NETWORK_CELLS", "induced,all"), ",")[[1]])
stopifnot(all(VARIANTS %in% c("induced", "all")))
variant_desc <- c(induced = "network built on INDUCED cells only", all = "network built on ALL cells (the original filter)")
set.seed(SEED)

source(UTILS_PATH)

# The C++ function is only a speed-up: R's own cor(method = "spearman") gives the same numbers (max difference ~2e-16)
# at ~2.5x the cost. If it cannot be compiled (e.g. a Mac without the gfortran 14.2 runtime CRAN's R links against),
# fall back to R. Force the fallback with COEXP_FORCE_R_CORRELATION=1.
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

# --------------------------------------------------------------------------
# data
# --------------------------------------------------------------------------
counts <- as.matrix(readRDS(file.path(PROCESSED_DIR, "scran_normalized_linear.rds")))
cat(sprintf("Counts: %d genes x %d cells (scran-normalized)\n", nrow(counts), ncol(counts)))

targetscan_all <- readRDS(TARGETS_PATH)
stopifnot(all(c("ensembl_gene_id", "Cumulative.weighted.context...score") %in% colnames(targetscan_all)))

# --------------------------------------------------------------------------
# fluorescence-to-cell join. IDENTICAL logic to fluorescence_accuracy_comparison.R: FACS index sort gives fluorescence
# per plate WELL; <prefix>_<plate>.well_barcodes.txt gives each cell's well. Validated empirically, because a wrong
# well geometry would scramble everything silently.
# --------------------------------------------------------------------------
norm_well <- function(x) {                      # "A01", "A1", "a1" -> "A1"; anything else -> NA
  x <- toupper(trimws(as.character(x))); ok <- grepl("^[A-Z]+[0-9]+$", x); out <- rep(NA_character_, length(x))
  out[ok] <- paste0(sub("^([A-Z]+)[0-9]+$", "\\1", x[ok]), as.integer(sub("^[A-Z]+([0-9]+)$", "\\1", x[ok]))); out
}
sp <- function(x, y) suppressWarnings(cor(x, y, method = "spearman", use = "complete.obs"))

build_fluorescence_table <- function() {
  meta <- readRDS(file.path(PROCESSED_DIR, "cell_metadata.rds"))
  need <- c("cell_id", "Dox_Group", "Dox_Concentration", "GFP_normalized")
  if (length(setdiff(need, colnames(meta))) > 0) stop("cell_metadata.rds lacks column(s): ", paste(setdiff(need, colnames(meta)), collapse = ", "))
  if (!"plate" %in% colnames(meta)) meta$plate <- sub("_.*", "", meta$cell_id)
  meta <- meta %>% transmute(cell_id, plate = as.character(plate), Dox_Group = as.character(Dox_Group),
                             Dox_Concentration = as.character(Dox_Concentration), GFP_normalized = as.numeric(GFP_normalized))
  cat(sprintf("cell_metadata: %d cells; Dox_Group: %s\n", nrow(meta), paste(names(table(meta$Dox_Group)), table(meta$Dox_Group), sep = "=", collapse = ", ")))

  read_fluor <- function(plate) {
    f <- file.path(RAW_DIR, sprintf("P%s_GFP_fluorescence.tsv", plate))
    if (!file.exists(f)) { cat(sprintf("*** plate %s: %s not found -- this plate is EXCLUDED ***\n", plate, basename(f))); return(NULL) }
    x <- read.delim(f, check.names = FALSE, stringsAsFactors = FALSE)
    if (!"GFP" %in% colnames(x)) stop(basename(f), " has no 'GFP' column; columns are: ", paste(colnames(x), collapse = ", "))
    if (!"Well_ID" %in% colnames(x)) {
      if (all(c("Row", "Column") %in% colnames(x))) x$Well_ID <- paste0(x$Row, x$Column)
      else stop(basename(f), " has neither 'Well_ID' nor 'Row'+'Column'; columns are: ", paste(colnames(x), collapse = ", "))
    }
    out <- tibble(plate = plate, well = norm_well(x$Well_ID), fluorescence = suppressWarnings(as.numeric(x$GFP)))
    out <- out[!is.na(out$well), ]
    n_dup <- sum(duplicated(out$well)); if (n_dup > 0) cat(sprintf("plate %s: %d wells appear more than once in the fluorescence file; keeping the first\n", plate, n_dup))
    out[!duplicated(out$well), ]
  }
  read_wellmap <- function(plate) {
    f <- file.path(RAW_DIR, sprintf("%s_%s.well_barcodes.txt", SEQ_PREFIX, plate))
    if (!file.exists(f)) stop("Well map not found: ", f, ". Files containing 'well_barcodes' in ", RAW_DIR, ": ", paste(list.files(RAW_DIR, pattern = "well_barcodes"), collapse = ", "))
    w <- read.delim(f, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
    if (!all(c("bc_set", "WellID") %in% colnames(w))) stop(basename(f), " must have columns bc_set and WellID; found: ", paste(colnames(w), collapse = ", "))
    tibble(cell_id = paste0(plate, "_", w$bc_set), well = norm_well(w$WellID))
  }
  fl <- bind_rows(lapply(PLATES, read_fluor))
  if (nrow(fl) == 0) stop("No fluorescence files found in ", RAW_DIR)
  wm <- bind_rows(lapply(unique(fl$plate), read_wellmap))
  cells <- meta %>% inner_join(wm, by = "cell_id") %>% left_join(fl, by = c("plate", "well"))

  cat("\n--- Join diagnostics per plate ---\n")
  dg <- cells %>% group_by(plate) %>% summarise(cells_with_well = n(), cells_with_fluorescence = sum(is.finite(fluorescence)), .groups = "drop") %>%
    mutate(match_frac = cells_with_fluorescence / cells_with_well)
  print(as.data.frame(dg), row.names = FALSE)
  bad <- dg$plate[dg$match_frac < MIN_MATCH_FRAC]
  if (length(bad) > 0) stop("Plate(s) ", paste(bad, collapse = ", "), ": fewer than ", 100 * MIN_MATCH_FRAC, "% of cells found a fluorescence value. Well IDs may be formatted differently.\n  example IDs in the well map: ",
                            paste(head(unique(wm$well), 5), collapse = ", "), "\n  example IDs in the fluorescence file: ", paste(head(unique(fl$well), 5), collapse = ", "))

  matched <- cells %>% dplyr::filter(is.finite(fluorescence), is.finite(GFP_normalized))
  rho_sanity <- sp(matched$fluorescence, matched$GFP_normalized)
  dose_levels <- c("Control", "0 Dox", "0.01 Dox", "0.1 Dox", "1 Dox")
  dose_med <- matched %>% dplyr::filter(Dox_Concentration %in% dose_levels) %>% group_by(Dox_Concentration) %>%
    summarise(n = n(), median_fluorescence = median(fluorescence), .groups = "drop") %>% arrange(match(Dox_Concentration, dose_levels))
  cat(sprintf("\n--- Join validation ---\nSpearman(fluorescence, GFP mRNA), all matched cells (n = %d): %.3f\n", nrow(matched), rho_sanity))
  cat("Median fluorescence by Dox concentration (expect a large rise from 0 Dox to the induced doses; a plateau between them is normal):\n")
  print(as.data.frame(dose_med), row.names = FALSE)
  if (!is.finite(rho_sanity) || rho_sanity < MIN_RHO_SANITY) {
    stop(sprintf("Fluorescence correlates only %.3f with GFP mRNA across matched cells (threshold %.1f); the well-to-cell join is almost certainly WRONG. ", rho_sanity, MIN_RHO_SANITY),
         "Do not use any result from this script until that is resolved.")
  }
  dmv <- setNames(dose_med$median_fluorescence, dose_med$Dox_Concentration)
  induced_doses <- intersect(c("0.01 Dox", "0.1 Dox", "1 Dox"), names(dmv))
  ok_rise <- "0 Dox" %in% names(dmv) && length(induced_doses) > 0 && dmv[[induced_doses[1]]] > max(2 * dmv[["0 Dox"]], dmv[["0 Dox"]] + 1e-9)
  step_ratio <- if (length(induced_doses) >= 2) dmv[induced_doses][-1] / dmv[induced_doses][-length(induced_doses)] else 1
  if (!ok_rise) cat("*** WARNING: median fluorescence of the lowest induced dose is not clearly (>2x) above 0 Dox -- check the join before trusting results. ***\n")
  else if (any(step_ratio < 0.75)) cat("*** WARNING: median fluorescence falls by more than 25% between induced doses -- check the join before trusting results. ***\n")
  else cat("PASSED: fluorescence tracks GFP mRNA and rises sharply with induction.\n")
  list(ind = matched %>% dplyr::filter(Dox_Group %in% INDUCED_GROUPS), sanity_rho = rho_sanity)
}

fl_tbl <- build_fluorescence_table()
ind <- fl_tbl$ind %>% dplyr::filter(cell_id %in% colnames(counts))
score_cells <- ind$cell_id
fluor_v <- ind$fluorescence
cat(sprintf("\nScoring cells: %d Dox-induced cells with a fluorescence value (no GFP-mRNA trimming: fluorescence is the ground truth and Spearman is rank-based)\n", length(score_cells)))
if (length(score_cells) < 50) stop("Too few induced cells with fluorescence (", length(score_cells), ").")

# --------------------------------------------------------------------------
# scoring functions (always: induced cells vs fluorescence)
# --------------------------------------------------------------------------
# Reference route: the established calculate_activity() pipeline.
rho_ref <- function(genes) {
  act <- calculate_activity(counts = counts, targets = genes)
  m <- inner_join(tibble(cell_id = names(act), activity = as.numeric(act)), tibble(cell_id = score_cells, fluor = fluor_v), by = "cell_id")
  m <- m[is.finite(m$activity) & is.finite(m$fluor), ]
  unname(cor.test(m$activity, m$fluor, method = "spearman", exact = FALSE)$estimate)
}
# Fast route for the thousands of null draws. activity = -log2(S/mean(S)) + 1 is strictly DEcreasing in the target sum S,
# so Spearman(activity, fluorescence) = -Spearman(S, fluorescence). Verified against rho_ref below before use.
rho_fast <- function(genes) {
  genes <- intersect(genes, rownames(counts))
  if (length(genes) < 2) return(NA_real_)
  s <- colSums(counts[genes, score_cells, drop = FALSE])
  ok <- is.finite(s) & s > 0
  if (sum(ok) < 10) return(NA_real_)
  -suppressWarnings(cor(s[ok], fluor_v[ok], method = "spearman"))
}

# the filter, as one reusable function
coexpr_select <- function(mat, thr = CORR_THRESHOLD, min_deg = MIN_DEGREE, keep_matrices = FALSE) {
  C <- spearman_full_cpp(mat); diag(C) <- 0
  A <- C > thr
  deg <- rowSums(A); wdeg <- rowSums(C * A)
  names(deg) <- names(wdeg) <- rownames(mat)
  out <- list(keep = names(deg)[deg >= min_deg], degree = deg, wdegree = wdeg, n_edges = sum(A) / 2)
  if (keep_matrices) { rownames(A) <- colnames(A) <- rownames(mat); out$adj <- A; out$cor <- C }
  out
}
top_k_by_degree <- function(sel, k) names(sel$degree)[order(-sel$degree, -sel$wdegree)][seq_len(k)]

# --------------------------------------------------------------------------
# [1] C++ vs R Spearman, once (all cells), at the level of the edge threshold
# --------------------------------------------------------------------------
cand_global <- unique(intersect(targetscan_all$ensembl_gene_id, rownames(counts)))
cand_global <- cand_global[apply(counts[cand_global, , drop = FALSE], 1, var) > 0]
cat("\n================ [1] C++ correlation validation ================\n")
val_ids <- sample(cand_global, min(300, length(cand_global)))   # drawn in both modes so the random stream is identical
if (USE_CPP) {
  Y <- rbind(c(1, 1, 1, 2, 3, 3), c(6, 5, 4, 3, 2, 1))
  ka_cpp <- spearman_full_cpp(Y)[1, 2]; ka_r <- cor(Y[1, ], Y[2, ], method = "spearman")
  cat(sprintf("known-answer case with ties: C++ %.10f vs R %.10f\n", ka_cpp, ka_r))
  Xv <- counts[val_ids, , drop = FALSE]
  Cc <- spearman_full_cpp(Xv); Cr <- cor(t(Xv), method = "spearman"); ut <- upper.tri(Cc)
  max_diff <- max(abs(Cc - Cr))
  n_discordant <- sum((Cc[ut] > CORR_THRESHOLD) != (Cr[ut] > CORR_THRESHOLD))
  cat(sprintf("%d real candidate genes: max |C++ - R| = %.2e; pairs on opposite sides of the %.1f threshold: %d of %d\n", length(val_ids), max_diff, CORR_THRESHOLD, n_discordant, sum(ut)))
  if (abs(ka_cpp - ka_r) > COR_TOL || max_diff > COR_TOL || n_discordant > 0) {
    stop("spearman_full_cpp() does not reproduce R's cor(method = 'spearman'). Make sure ", CPP_PATH, " is the CORRECTED fullSpearman.cpp (average-rank ties).")
  }
  cat("PASSED: the C++ function agrees with R's Spearman, including at the edge threshold.\n")
} else {
  max_diff <- NA_real_; n_discordant <- NA_integer_
  cat("C++ not in use: all correlations are computed with R's cor(method = 'spearman'), the reference implementation, so there is nothing to validate here.\n")
}

# scoring cross-check, once
ts_sorted <- targetscan_all[order(targetscan_all$Cumulative.weighted.context...score), ]
top200_genes <- intersect(ts_sorted[seq_len(N_TOP_TARGETS), ]$ensembl_gene_id, rownames(counts))
ts_score_gene <- tapply(targetscan_all$Cumulative.weighted.context...score, targetscan_all$ensembl_gene_id, min, na.rm = TRUE)
chk_sets <- list(top200 = top200_genes, all_candidates = cand_global, random_a = sample(cand_global, 35), random_b = sample(cand_global, 35))
chk <- sapply(chk_sets, function(s) c(ref = rho_ref(s), fast = rho_fast(s)))
cat("\nScoring cross-check (established calculate_activity() pipeline vs the fast route used for the null draws), induced cells vs fluorescence:\n"); print(round(t(chk), 6))
if (max(abs(chk["ref", ] - chk["fast", ])) > 1e-8) stop("Fast scoring route does not reproduce calculate_activity()-based rho; do not use the null results.")
cat("PASSED: both routes agree.\n")

# ===========================================================================
# one network variant: observed filter, [2] random subsets, [3] random pools, [4] features
# ===========================================================================
run_variant <- function(variant) {
  cat(sprintf("\n################################################################\n# VARIANT: %s\n################################################################\n", variant_desc[[variant]]))
  net_cells <- if (variant == "induced") score_cells else colnames(counts)
  cn <- counts[, net_cells, drop = FALSE]
  cat(sprintf("Network cells: %d | scoring cells (induced, with fluorescence): %d\n", length(net_cells), length(score_cells)))

  cand_all <- unique(intersect(targetscan_all$ensembl_gene_id, rownames(cn)))
  cv <- apply(cn[cand_all, , drop = FALSE], 1, var)
  cand_ids <- cand_all[is.finite(cv) & cv > 0]
  cat(sprintf("TargetScan candidates present in the scran matrix: %d (non-zero variance in these cells: %d); exact-zero fraction %.3f\n",
              length(cand_all), length(cand_ids), mean(cn[cand_ids, ] == 0)))

  # ---------------- observed ----------------
  cand_mat <- cn[cand_ids, , drop = FALSE]
  obs <- coexpr_select(cand_mat, keep_matrices = TRUE)
  K <- length(obs$keep)
  thr_edges <- sapply(c(0.2, 0.3, 0.4, 0.5), function(t) sum(obs$cor > t) / 2)
  cat(sprintf("Edges at rho > 0.2 / 0.3 / 0.4 / 0.5: %s\n", paste(thr_edges, collapse = " / ")))
  cat(sprintf("Candidates: %d | edges (rho > %.2f): %d | genes with degree >= %d: %d | degree min=%d, median=%.0f, max=%d\n",
              length(cand_ids), CORR_THRESHOLD, obs$n_edges, MIN_DEGREE, K, min(obs$degree), median(obs$degree), max(obs$degree)))
  if (K < MIN_SURV) {
    cat(sprintf("*** Only %d gene(s) survive (< %d) in this variant; nothing further to validate. Within these cells co-expression among the candidates is weak at this threshold -- see the edge counts above for other thresholds (COEXP_CORR_THRESHOLD). ***\n", K, MIN_SURV))
    return(list(variant = variant, skipped = TRUE, K = K, n_edges = obs$n_edges, thr_edges = thr_edges, selected = obs$keep, n_net_cells = length(net_cells)))
  }
  g <- igraph::graph_from_adjacency_matrix(1 * obs$adj, mode = "undirected", diag = FALSE)
  stopifnot(all(igraph::degree(g)[names(obs$degree)] == obs$degree))

  rho_obs     <- rho_fast(obs$keep)
  rho_all     <- rho_fast(cand_ids)
  rho_top200  <- rho_fast(top200_genes)
  ts_topK     <- names(sort(ts_score_gene[cand_ids]))[seq_len(K)]      # TargetScan's own top-K at the SAME size
  rho_ts_topK <- rho_fast(ts_topK)
  cat(sprintf("\nrho vs GFP FLUORESCENCE (induced cells):\n  co-expression-selected (K = %d genes)      %.3f\n  TargetScan top-%d, same size                 %.3f\n  TargetScan top 200 (%d genes present)        %.3f\n  all %d candidates                          %.3f\n",
              K, rho_obs, K, rho_ts_topK, length(top200_genes), rho_top200, length(cand_ids), rho_all))

  # ---------------- [2] random subsets ----------------
  cat(sprintf("\n================ [2] Random subsets of the candidate pool (K = %d, n = %d) ================\n", K, N_NULL))
  null_rho <- replicate(N_NULL, rho_fast(sample(cand_ids, K)))
  p_null <- (1 + sum(null_rho >= rho_obs, na.rm = TRUE)) / (1 + sum(!is.na(null_rho)))
  cat(sprintf("random K-subsets: median %.3f | 95th pct %.3f | 99th pct %.3f | max %.3f\n",
              median(null_rho, na.rm = TRUE), quantile(null_rho, 0.95, na.rm = TRUE), quantile(null_rho, 0.99, na.rm = TRUE), max(null_rho, na.rm = TRUE)))
  cat(sprintf("observed %.3f exceeds %.1f%% of random subsets; empirical p = %.4f (one-sided; smallest attainable %.4f)\n", rho_obs, 100 * mean(null_rho < rho_obs, na.rm = TRUE), p_null, 1 / (1 + N_NULL)))
  cat(sprintf("Equal-size comparison with TargetScan's own top-%d: co-expression %.3f vs TargetScan %.3f (difference %+.3f); the two sets share %d of %d genes\n",
              K, rho_obs, rho_ts_topK, rho_obs - rho_ts_topK, length(intersect(obs$keep, ts_topK)), K))

  # ---------------- [3] random pools ----------------
  cat(sprintf("\n================ [3] Random-pool control (%d expression-matched pools) ================\n", N_POOLS))
  var_n <- apply(cn, 1, var); mean_n <- rowMeans(cn); det_n <- rowMeans(cn > 0)      # computed on the NETWORK cells
  bin_n <- cut(rank(mean_n, ties.method = "average") / length(mean_n), breaks = seq(0, 1, length.out = N_EXPR_BINS + 1), include.lowest = TRUE, labels = FALSE)
  names(bin_n) <- rownames(cn)
  nontarget_ids <- setdiff(rownames(cn), targetscan_all$ensembl_gene_id)               # excludes ALL TargetScan candidates
  nontarget_ids <- nontarget_ids[is.finite(var_n[nontarget_ids]) & var_n[nontarget_ids] > 0]
  need <- table(factor(bin_n[cand_ids], levels = seq_len(N_EXPR_BINS)))
  pool_by_bin <- split(nontarget_ids, factor(bin_n[nontarget_ids], levels = seq_len(N_EXPR_BINS)))
  shortfall <- sum(pmax(as.integer(need) - lengths(pool_by_bin), 0))
  cat(sprintf("Non-target gene pool: %d genes; each random pool mirrors the candidates' expression distribution (%d bins, matched on the network cells).\n", length(nontarget_ids), N_EXPR_BINS))
  if (shortfall > 0) cat(sprintf("NOTE: %d candidate slots have no available non-target gene in their bin; those pools will be slightly smaller.\n", shortfall))
  draw_matched_pool <- function() unlist(lapply(seq_len(N_EXPR_BINS), function(b) {
    n_need <- as.integer(need[b]); avail <- pool_by_bin[[b]]
    if (n_need == 0 || length(avail) == 0) return(NULL)
    sample(avail, min(n_need, length(avail)))
  }), use.names = FALSE)

  rho_obs_topK <- rho_fast(top_k_by_degree(obs, K))
  pool_res <- vector("list", N_POOLS)
  for (b in seq_len(N_POOLS)) {
    pool <- draw_matched_pool()
    sel <- coexpr_select(cn[pool, , drop = FALSE])
    # plain variables first, so no tibble column can shadow the gene vector `pool`
    pool_n_genes <- length(pool); pool_med_expr <- median(mean_n[pool]); pool_med_det <- median(det_n[pool])
    pool_rho_surv <- if (length(sel$keep) >= MIN_SURV) rho_fast(sel$keep) else NA_real_
    pool_rho_topK <- rho_fast(top_k_by_degree(sel, K))
    pool_res[[b]] <- tibble(pool_id = b, n_genes = pool_n_genes, n_edges = sel$n_edges, n_survivors = length(sel$keep),
                            rho_survivors = pool_rho_surv, rho_topK = pool_rho_topK, med_expr = pool_med_expr, med_det = pool_med_det)
    if (b %% 20 == 0) cat(sprintf("  ... %d/%d pools done\n", b, N_POOLS))
  }
  pool_res <- bind_rows(pool_res)
  cat(sprintf("\nBalance check (candidates vs matched pools): median mean-expression %.3f vs %.3f | median fraction detected %.3f vs %.3f | pool size %d vs %d candidates\n",
              median(mean_n[cand_ids]), mean(pool_res$med_expr), median(det_n[cand_ids]), mean(pool_res$med_det), round(mean(pool_res$n_genes)), length(cand_ids)))
  bal_rel <- abs(mean(pool_res$med_expr) - median(mean_n[cand_ids])) / median(mean_n[cand_ids])
  if (bal_rel > 0.15) cat(sprintf("*** WARNING: matched pools differ from the candidates' median expression by %.0f%% -- the control is not well matched; do not interpret it. ***\n", 100 * bal_rel))
  cat(sprintf("Edges (rho > %.2f):   observed candidates %d | random pools median %.0f (range %d-%d)\n", CORR_THRESHOLD, obs$n_edges, median(pool_res$n_edges), min(pool_res$n_edges), max(pool_res$n_edges)))
  cat(sprintf("Survivors (degree >= %d): observed %d | random pools median %.0f (range %d-%d); %d of %d pools have >= %d\n",
              MIN_DEGREE, K, median(pool_res$n_survivors), min(pool_res$n_survivors), max(pool_res$n_survivors), sum(pool_res$n_survivors >= MIN_SURV), N_POOLS, MIN_SURV))
  cat(sprintf("rho of the top-%d genes by degree (size-matched), vs fluorescence: observed %.3f | random pools median %.3f, range %.3f to %.3f\n",
              K, rho_obs_topK, median(pool_res$rho_topK, na.rm = TRUE), min(pool_res$rho_topK, na.rm = TRUE), max(pool_res$rho_topK, na.rm = TRUE)))
  cat("(the sign matters: strongly co-expressed random genes can correlate with the score in EITHER direction; a large value of either sign means the filter finds generic structure.\n",
      " Fluorescence also scales with cell size, so a size-linked axis can appear here exactly as it did for GFP mRNA.)\n", sep = "")
  p_pool_rho   <- (1 + sum(pool_res$rho_topK >= rho_obs_topK, na.rm = TRUE)) / (1 + sum(!is.na(pool_res$rho_topK)))
  p_pool_edges <- (1 + sum(pool_res$n_edges >= obs$n_edges)) / (1 + N_POOLS)
  cat(sprintf("empirical p: rho of top-%d by degree >= observed in random pools: %.4f | edges >= observed: %.4f   (smallest attainable with %d pools: %.4f)\n",
              K, p_pool_rho, p_pool_edges, N_POOLS, 1 / (1 + N_POOLS)))

  # ---------------- [4] selected vs other candidates ----------------
  cat("\n================ [4] Selected genes vs other candidates ================\n")
  feat <- tibble(ensembl_gene_id = cand_ids, mean_expr = rowMeans(cand_mat), detection = rowMeans(cand_mat > 0), selected = cand_ids %in% obs$keep)
  feat$context_score <- as.numeric(ts_score_gene[cand_ids])
  mer8_cols <- grep("8mer", colnames(targetscan_all), value = TRUE)
  if (length(mer8_cols) > 0) {
    tmp <- targetscan_all[targetscan_all$ensembl_gene_id %in% cand_ids, ]
    M <- do.call(cbind, lapply(mer8_cols, function(cn_) suppressWarnings(as.numeric(tmp[[cn_]]))))
    tmp$mer8 <- rowSums(M, na.rm = TRUE)
    feat$mer8_sites <- as.numeric(tapply(tmp$mer8, tmp$ensembl_gene_id, max)[cand_ids])
  } else cat("No column containing '8mer' found in the TargetScan table; 8mer comparison skipped.\n")
  if (file.exists(UTR_PATH)) {
    u <- readRDS(UTR_PATH)
    if (all(c("gid", "nchar") %in% colnames(u))) feat$utr_length <- as.numeric(tapply(u$nchar, u$gid, max)[cand_ids])
    else cat("UTR file lacks gid/nchar columns; UTR comparison skipped.\n")
  } else cat("UTR file not found; UTR comparison skipped.\n")
  compare_feature <- function(col, label) {
    if (!col %in% colnames(feat)) return(NULL)
    x <- feat[[col]][feat$selected]; y <- feat[[col]][!feat$selected]
    if (sum(is.finite(x)) < 3 || sum(is.finite(y)) < 3) return(NULL)
    tibble(feature = label, median_selected = median(x, na.rm = TRUE), median_other = median(y, na.rm = TRUE), wilcoxon_p = suppressWarnings(wilcox.test(x, y)$p.value))
  }
  feat_cmp <- bind_rows(compare_feature("context_score", "TargetScan context score (more negative = stronger)"), compare_feature("mer8_sites", "8mer sites"),
                        compare_feature("utr_length", "3'UTR length (nt)"), compare_feature("mean_expr", "mean expression (scran)"), compare_feature("detection", "fraction of cells detected"))
  cat("Median per group (selected vs the other candidates), two-sided Wilcoxon:\n"); print(as.data.frame(feat_cmp), digits = 4, row.names = FALSE)
  feat$ts_rank <- rank(feat$context_score, ties.method = "first")
  enrich <- bind_rows(lapply(c(100, 200, 500), function(X) {
    top <- feat$ts_rank <= X; k <- sum(feat$selected); m <- sum(top); x <- sum(top & feat$selected)
    tibble(top_X_by_context_score = X, selected_in_top = x, expected_by_chance = k * m / nrow(feat), hypergeometric_p = phyper(x - 1, m, nrow(feat) - m, k, lower.tail = FALSE))
  }))
  cat(sprintf("\nOverlap of the %d selected genes with TargetScan's strongest candidates (ranked among the %d present):\n", K, nrow(feat))); print(as.data.frame(enrich), digits = 4, row.names = FALSE)

  # ---------------- plots ----------------
  sub_txt <- sprintf("%s; scored against GFP fluorescence on %d induced cells", variant_desc[[variant]], length(score_cells))
  refs <- tibble(label = c(sprintf("co-expression selected (K=%d)", K), sprintf("TargetScan top-%d (same size)", K), sprintf("TargetScan top 200 (%d genes)", length(top200_genes)), sprintf("all %d candidates", length(cand_ids))),
                 value = c(rho_obs, rho_ts_topK, rho_top200, rho_all))
  p1 <- ggplot(tibble(rho = null_rho[!is.na(null_rho)]), aes(rho)) + geom_histogram(bins = 40, fill = "grey75", color = "grey40") +
    geom_vline(data = refs, aes(xintercept = value, color = label, linetype = label), linewidth = 1.1) +
    theme_bw(base_size = 13) + theme(legend.position = "bottom", legend.direction = "vertical", legend.title = element_blank(), plot.subtitle = element_text(size = 9)) +
    labs(x = "Spearman rho vs GFP fluorescence (induced cells)", y = "random subsets", title = sprintf("Random %d-gene subsets of the candidate pool (n = %d)", K, N_NULL),
         subtitle = sprintf("%s\nempirical p for the co-expression-selected set = %.4f", sub_txt, p_null))
  ggsave(file.path(OUT_DIR, sprintf("coexp_fluor_%s_random_subsets.pdf", variant)), p1, width = 8, height = 7.5)

  pool_long <- pool_res %>% transmute(`rho of top-K genes by degree` = rho_topK, `edges above threshold` = n_edges, `survivors (degree >= min)` = n_survivors) %>%
    pivot_longer(everything(), names_to = "metric", values_to = "value")
  obs_lines <- tibble(metric = c("rho of top-K genes by degree", "edges above threshold", "survivors (degree >= min)"), value = c(rho_obs_topK, obs$n_edges, K))
  p2 <- ggplot(pool_long, aes(value)) + geom_histogram(bins = 25, fill = "grey75", color = "grey40", na.rm = TRUE) +
    geom_vline(data = obs_lines, aes(xintercept = value), color = "red", linewidth = 1.1) + facet_wrap(~metric, scales = "free", ncol = 1) +
    theme_bw(base_size = 13) + theme(plot.subtitle = element_text(size = 9)) +
    labs(x = NULL, y = "random pools", title = "Expression-matched random pools", subtitle = sprintf("n = %d pools; red line = observed candidate pool\n%s", N_POOLS, sub_txt))
  ggsave(file.path(OUT_DIR, sprintf("coexp_fluor_%s_random_pools.pdf", variant)), p2, width = 7, height = 9)

  feat_long <- feat %>% transmute(selected = ifelse(selected, "selected", "other candidates"), `context score` = context_score,
      `log10 mean expression` = log10(mean_expr + 1e-6), `fraction detected` = detection,
      `log10 3'UTR length` = if ("utr_length" %in% colnames(feat)) log10(utr_length) else NA_real_,
      `8mer sites` = if ("mer8_sites" %in% colnames(feat)) mer8_sites else NA_real_) %>%
    pivot_longer(-selected, names_to = "feature", values_to = "value") %>% dplyr::filter(is.finite(value))
  p3 <- ggplot(feat_long, aes(selected, value, fill = selected)) + geom_boxplot(outlier.alpha = 0.3, show.legend = FALSE) + facet_wrap(~feature, scales = "free_y") +
    theme_bw(base_size = 12) + theme(plot.subtitle = element_text(size = 9)) + labs(x = NULL, y = NULL, title = "Selected genes vs other candidates", subtitle = variant_desc[[variant]])
  ggsave(file.path(OUT_DIR, sprintf("coexp_fluor_%s_features.pdf", variant)), p3, width = 9, height = 6)

  list(variant = variant, skipped = FALSE, n_net_cells = length(net_cells), K = K, selected = obs$keep, degree = obs$degree, n_edges = obs$n_edges, thr_edges = thr_edges,
       rho = c(selected = rho_obs, ts_topK = rho_ts_topK, top200 = rho_top200, all_candidates = rho_all),
       null_rho = null_rho, p_null = p_null, pool_res = pool_res, rho_obs_topK = rho_obs_topK, p_pool_rho = p_pool_rho, p_pool_edges = p_pool_edges,
       features = feat, feature_comparison = feat_cmp, enrichment = enrich)
}

results <- lapply(VARIANTS, run_variant)
names(results) <- VARIANTS

# --------------------------------------------------------------------------
# cross-variant summary
# --------------------------------------------------------------------------
cat("\n################################################################\n# SUMMARY ACROSS NETWORK VARIANTS (all scored vs fluorescence, induced cells)\n################################################################\n")
summ <- bind_rows(lapply(results, function(r) {
  if (isTRUE(r$skipped)) return(tibble(network = r$variant, network_cells = r$n_net_cells, K = r$K, rho_selected = NA_real_, rho_TargetScan_topK = NA_real_, p_random_subsets = NA_real_,
                                       pools_median_rho = NA_real_, p_pools_rho = NA_real_, p_pools_edges = NA_real_))
  tibble(network = r$variant, network_cells = r$n_net_cells, K = r$K, rho_selected = r$rho[["selected"]], rho_TargetScan_topK = r$rho[["ts_topK"]],
         p_random_subsets = r$p_null, pools_median_rho = median(r$pool_res$rho_topK, na.rm = TRUE), p_pools_rho = r$p_pool_rho, p_pools_edges = r$p_pool_edges)
}))
print(as.data.frame(summ %>% mutate(across(c(starts_with("rho"), ends_with("median_rho")), ~ round(.x, 3)), across(starts_with("p_"), ~ signif(.x, 2)))), row.names = FALSE)
if (all(c("induced", "all") %in% names(results))) {
  ov <- length(intersect(results$induced$selected, results$all$selected))
  cat(sprintf("\nSelected sets: induced-network K = %d, all-cell-network K = %d; they share %d genes\n", results$induced$K, results$all$K, ov))
  cat("(a small overlap means the induction step / Dox exposure, not within-induced co-variation, shaped the all-cell selection)\n")
}

saveRDS(list(variants = results, summary = summ,
             scoring = list(ground_truth = "GFP fluorescence (FITC-A) from FACS index sort", cells = "Dox-induced cells", n_cells = length(score_cells),
                            induced_groups = INDUCED_GROUPS, fluorescence_vs_mRNA_all_cells = fl_tbl$sanity_rho),
             params = list(corr_threshold = CORR_THRESHOLD, min_degree = MIN_DEGREE, n_null = N_NULL, n_pools = N_POOLS, used_cpp = USE_CPP,
                           cpp_validation = list(max_diff = max_diff, discordant_pairs = n_discordant))),
        file.path(OUT_DIR, "res_coexpression_validation_fluorescence.rds"))

cat("\n================ How to read this ================\n")
cat("Two network variants are scored identically (induced cells vs fluorescence), so any difference between them is due to which genes were selected.\n",
    "[2] If the observed rho sits inside the bulk of the random-subset histogram, selection adds nothing beyond picking candidates at random. Also compare\n",
    "    with TargetScan's own top-K at the same size.\n",
    "[3] If expression-matched random pools often give top-K genes with rho comparable to the observed one, the filter mostly finds generic co-expression\n",
    "    (e.g. a cell-size axis -- fluorescence scales with cell size too), not miR-124 targets. LIMIT: this control cannot tell miRNA-driven co-expression\n",
    "    apart from co-expression due to other shared properties of TargetScan candidates (3'UTR length, tissue-specific modules).\n",
    "[4] Selected genes that look like ordinary strong candidates support a target interpretation; mainly highly expressed or long-UTR genes point to a\n",
    "    technical or compositional explanation.\n",
    "'induced' vs 'all': if the all-cell network scores clearly lower than the induced-cell network under induced-cell scoring, its selection was driven\n",
    "    partly by the induction step (or Dox exposure), which carries no within-induced signal.\n", sep = "")
cat(sprintf("\nSaved plots (coexp_fluor_<variant>_*.pdf) and res_coexpression_validation_fluorescence.rds to %s\n", OUT_DIR))

