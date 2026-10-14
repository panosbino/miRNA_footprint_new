suppressMessages({
  library(dplyr); library(tibble); library(tidyr); library(ggplot2)
  library(Rcpp); library(RcppArmadillo)
})

# ===========================================================================
# How do the co-expressiondplyr::filter's two cutoffs relate to its accuracy?
# Accuracy = Spearman rho between the target-sum activity score of thedplyr::selected
# genes and GFP FLUORESCENCE, on Dox-induced cells only.
#
# Thedplyr::filter has two parameters. Each is swept in two ways:
#   rho cutoff        : (a) specific numbers      edge if rho > t        t = 0.20 ... 0.60
#                       (b) percentages           edge if the pair is in the top p% of all
#                                                 candidate-pair correlations (adapts to the
#                                                 network, so induced and all-cell networks
#                                                 are comparable)
#   edges per gene   : (a) specific numbers      keep genes with >= d edges
#                       (b) percentages           keep the top q% of genes by connectivity
#                                                 (fixes the number of genes::selected)
# All four combinations are evaluated. For every grid point the table reports the
# number of edges/genes both as counts and as percentages, plus three yardsticks that
# a bare accuracy cannot give:
#   * TargetScan's own top-K at the same K (does the filter beat the ranking it started from?)
#   * random K-subsets of the candidates
#   * the IDENTICAL rule applied to expression-matched NON-target pools (what the
#   filter reaches on generic co-expression structure alone; same K as the candidates)
#
# Two network variants, as in coexpression_filter_validation_fluorescence.R:
#   "induced": network from induced cells only   "all": network from all cells.
# Scoring is always induced cells vs fluorescence. Fluorescence is never used to select genes.
#
# CAUTION: choosing the best cutoff by this table and then quoting its accuracy would be
# optimistic (tuned on the same cells it is scored on). Read the table as a sensitivity
# analysis; a cutoff to be reported should be chosen on separate data (e.g. one plate held out).
# ===========================================================================

env_int <- function(name, default) as.integer(Sys.getenv(name, as.character(default)))

BASE_DIR      <- Sys.getenv("COEXP_BASE_DIR", Sys.getenv("HEK_BASE_DIR", "~/Desktop/Projects/miRNA_footprint_new"))
RAW_DIR       <- file.path(BASE_DIR, "datasets/HEK_SS3/raw")
PROCESSED_DIR <- file.path(BASE_DIR, "datasets/HEK_SS3/processed")
OUT_DIR       <- file.path(BASE_DIR, "analysis/HEK_SS3/comparisons")
TARGETS_PATH  <- file.path(BASE_DIR, "resources/human/Targets__combined_124-3p_124-3p.2_506-3p.rds")
UTILS_PATH    <- file.path(BASE_DIR, "scripts/Utils.R")          # defines calculate_activity()
CPP_PATH      <- file.path(BASE_DIR, "scripts/fullSpearman.cpp") # MUST be the corrected version
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

SEQ_PREFIX     <- "P32156"
PLATES         <- c("101", "103", "105")
INDUCED_GROUPS <- "Dox Induced"
MIN_MATCH_FRAC <- 0.5
MIN_RHO_SANITY <- 0.2

# ---- the grid (edit freely) ----
RHO_CUTOFFS <- c(0.20, 0.25, 0.30, 0.35, 0.40, 0.45, 0.50, 0.55, 0.60)   # specific rho cutoffs
MIN_DEGREES <- c(1, 2, 3, 4, 5, 7, 10, 15, 20)                            # specific minimum numbers of edges per gene
EDGE_PCTS   <- c(0.005, 0.01, 0.02, 0.05, 0.1, 0.25, 0.5, 1)              # % of all candidate pairs kept as edges
DEG_PCTS    <- c(1, 2, 3, 5, 10, 20, 30)                                  # % of candidates kept (most connected first)

MIN_K        <- 3                                    # grid points::selecting fewer genes are reported as NA
N_POOLS_GRID <- env_int("COEXP_N_POOLS_GRID", 50)    # random pools run through EVERY grid point (smallest attainable p = 1/(n+1))
N_NULL_GRID  <- env_int("COEXP_N_NULL_GRID", 200)    # random K-subsets per distinct K
N_TOP_TARGETS <- 200
N_EXPR_BINS  <- 30
COR_TOL      <- 1e-6
SEED         <- 2026
VARIANTS     <- trimws(strsplit(Sys.getenv("COEXP_NETWORK_CELLS", "induced,all"), ",")[[1]])
stopifnot(all(VARIANTS %in% c("induced", "all")))
variant_desc <- c(induced = "network built on INDUCED cells only", all = "network built on ALL cells (the original::filter)")
set.seed(SEED)
options(width = 200)

source(UTILS_PATH)

USE_CPP <- FALSE
if (Sys.getenv("COEXP_FORCE_R_CORRELATION", "0") != "1") {
  USE_CPP <- tryCatch({ sourceCpp(CPP_PATH); TRUE }, error = function(e) {
    cat("\n*** Could not compile/link ", CPP_PATH, " (", conditionMessage(e), ").\n",
        "    Falling back to R's cor(method = 'spearman'): identical results, a bit slower. ***\n\n", sep = "")
    FALSE
  })
}
if (!USE_CPP) spearman_full_cpp <- function(X) cor(t(X), method = "spearman")   # genes x cells in, genes x genes out

# --------------------------------------------------------------------------
# data and the fluorescence join (IDENTICAL logic to the other fluorescence scripts)
# --------------------------------------------------------------------------
counts <- as.matrix(readRDS(file.path(PROCESSED_DIR, "scran_normalized_linear.rds")))
cat(sprintf("Counts: %d genes x %d cells (scran-normalized)\n", nrow(counts), ncol(counts)))
targetscan_all <- readRDS(TARGETS_PATH)
stopifnot(all(c("ensembl_gene_id", "Cumulative.weighted.context...score") %in% colnames(targetscan_all)))

norm_well <- function(x) {
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
  read_fluor <- function(plate) {
    f <- file.path(RAW_DIR, sprintf("P%s_GFP_fluorescence.tsv", plate))
    if (!file.exists(f)) { cat(sprintf("*** plate %s: %s not found -- this plate is EXCLUDED ***\n", plate, basename(f))); return(NULL) }
    x <- read.delim(f, check.names = FALSE, stringsAsFactors = FALSE)
    if (!"GFP" %in% colnames(x)) stop(basename(f), " has no 'GFP' column; columns are: ", paste(colnames(x), collapse = ", "))
    if (!"Well_ID" %in% colnames(x)) {
      if (all(c("Row", "Column") %in% colnames(x))) x$Well_ID <- paste0(x$Row, x$Column)
      else stop(basename(f), " has neither 'Well_ID' nor 'Row'+'Column'")
    }
    out <- tibble(plate = plate, well = norm_well(x$Well_ID), fluorescence = suppressWarnings(as.numeric(x$GFP)))
    out <- out[!is.na(out$well), ]; out[!duplicated(out$well), ]
  }
  read_wellmap <- function(plate) {
    f <- file.path(RAW_DIR, sprintf("%s_%s.well_barcodes.txt", SEQ_PREFIX, plate))
    if (!file.exists(f)) stop("Well map not found: ", f)
    w <- read.delim(f, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
    if (!all(c("bc_set", "WellID") %in% colnames(w))) stop(basename(f), " must have columns bc_set and WellID")
    tibble(cell_id = paste0(plate, "_", w$bc_set), well = norm_well(w$WellID))
  }
  fl <- bind_rows(lapply(PLATES, read_fluor)); if (nrow(fl) == 0) stop("No fluorescence files found in ", RAW_DIR)
  wm <- bind_rows(lapply(unique(fl$plate), read_wellmap))
  cells <- meta %>% inner_join(wm, by = "cell_id") %>% left_join(fl, by = c("plate", "well"))
  dg <- cells %>% group_by(plate) %>% summarise(n = n(), with_fluor = sum(is.finite(fluorescence)), .groups = "drop") %>% mutate(frac = with_fluor / n)
  if (any(dg$frac < MIN_MATCH_FRAC)) stop("Plate(s) ", paste(dg$plate[dg$frac < MIN_MATCH_FRAC], collapse = ", "), ": too few cells found a fluorescence value; well IDs may be formatted differently.")
  matched <- cells %>%dplyr::filter(is.finite(fluorescence), is.finite(GFP_normalized))
  rho_sanity <- sp(matched$fluorescence, matched$GFP_normalized)
  cat(sprintf("Join check: %d cells matched; Spearman(fluorescence, GFP mRNA) = %.3f\n", nrow(matched), rho_sanity))
  if (!is.finite(rho_sanity) || rho_sanity < MIN_RHO_SANITY) stop(sprintf("Fluorescence correlates only %.3f with GFP mRNA; the well-to-cell join is almost certainly WRONG. Do not use any result from this script.", rho_sanity))
  matched %>%dplyr::filter(Dox_Group %in% INDUCED_GROUPS)
}
ind <- build_fluorescence_table() %>%dplyr::filter(cell_id %in% colnames(counts))
score_cells <- ind$cell_id; fluor_v <- ind$fluorescence
cat(sprintf("Scoring cells: %d Dox-induced cells with a fluorescence value\n", length(score_cells)))
if (length(score_cells) < 50) stop("Too few induced cells with fluorescence (", length(score_cells), ").")

# --------------------------------------------------------------------------
# scoring: rho of the target-sum activity vs fluorescence on induced cells
# --------------------------------------------------------------------------
rho_ref <- function(genes) {                      # established calculate_activity() route
  act <- calculate_activity(counts = counts, targets = genes)
  m <- inner_join(tibble(cell_id = names(act), activity = as.numeric(act)), tibble(cell_id = score_cells, fluor = fluor_v), by = "cell_id")
  m <- m[is.finite(m$activity) & is.finite(m$fluor), ]
  unname(cor.test(m$activity, m$fluor, method = "spearman", exact = FALSE)$estimate)
}
rho_fast <- function(genes) {                     # activity is strictly decreasing in the target sum S, so rho = -Spearman(S, fluorescence)
  genes <- intersect(genes, rownames(counts))
  if (length(genes) < 2) return(NA_real_)
  s <- colSums(counts[genes, score_cells, drop = FALSE]); ok <- is.finite(s) & s > 0
  if (sum(ok) < 10) return(NA_real_)
  -suppressWarnings(cor(s[ok], fluor_v[ok], method = "spearman"))
}

# --------------------------------------------------------------------------
# checks, once: C++ vs R at the edge-threshold level; fast scoring vs calculate_activity()
# --------------------------------------------------------------------------
cand_global <- unique(intersect(targetscan_all$ensembl_gene_id, rownames(counts)))
cand_global <- cand_global[apply(counts[cand_global, , drop = FALSE], 1, var) > 0]
val_ids <- sample(cand_global, min(300, length(cand_global)))
if (USE_CPP) {
  Y <- rbind(c(1, 1, 1, 2, 3, 3), c(6, 5, 4, 3, 2, 1))
  Xv <- counts[val_ids, , drop = FALSE]; Cc <- spearman_full_cpp(Xv); Cr <- cor(t(Xv), method = "spearman"); ut <- upper.tri(Cc)
  ok <- abs(spearman_full_cpp(Y)[1, 2] - cor(Y[1, ], Y[2, ], method = "spearman")) <= COR_TOL && max(abs(Cc - Cr)) <= COR_TOL &&
        sum((Cc[ut] > 0.4) != (Cr[ut] > 0.4)) == 0
  if (!ok) stop("spearman_full_cpp() does not reproduce R's cor(method = 'spearman'); make sure ", CPP_PATH, " is the CORRECTED fullSpearman.cpp.")
  cat("C++ correlation check PASSED.\n")
} else cat("C++ not in use: correlations computed with R's cor(method = 'spearman').\n")
ts_sorted <- targetscan_all[order(targetscan_all$Cumulative.weighted.context...score), ]
top200_genes <- intersect(ts_sorted[seq_len(N_TOP_TARGETS), ]$ensembl_gene_id, rownames(counts))
ts_score_gene <- tapply(targetscan_all$Cumulative.weighted.context...score, targetscan_all$ensembl_gene_id, min, na.rm = TRUE)
chk <- sapply(list(top200 = top200_genes, all_candidates = cand_global, random_a = sample(cand_global, 35)), function(s) c(ref = rho_ref(s), fast = rho_fast(s)))
if (max(abs(chk["ref", ] - chk["fast", ])) > 1e-8) stop("Fast scoring route does not reproduce calculate_activity()-based rho.")
cat("Scoring cross-check PASSED.\n")

# --------------------------------------------------------------------------
# the grid and the network machinery
# --------------------------------------------------------------------------
edge_rules <- bind_rows(tibble(edge_type = "rho", edge_value = RHO_CUTOFFS), tibble(edge_type = "pct", edge_value = EDGE_PCTS)) %>% mutate(edge_id = paste0(edge_type, ":", edge_value))
deg_rules  <- bind_rows(tibble(deg_type = "min", deg_value = MIN_DEGREES), tibble(deg_type = "pct", deg_value = DEG_PCTS)) %>% mutate(deg_id = paste0(deg_type, ":", deg_value))
gdf <- as_tibble(merge(edge_rules, deg_rules, by = NULL)) %>% mutate(grid_id = row_number())
cat(sprintf("Grid: %d edge rules x %d degree rules = %d settings per network variant\n", nrow(edge_rules), nrow(deg_rules), nrow(gdf)))

# Everything about one correlation matrix C (diagonal 0, rownames = genes), for ALL edge rules at once:
# edges, degree, a ranking of genes by connectivity (ties broken at random, never by input order).
network_stats <- function(C) {
  n <- nrow(C); npairs <- n * (n - 1) / 2
  vals <- sort(C[upper.tri(C)], decreasing = TRUE)
  thr <- c(setNames(RHO_CUTOFFS, paste0("rho:", RHO_CUTOFFS)),
           setNames(vals[pmin(npairs, pmax(1, ceiling(EDGE_PCTS / 100 * npairs)))], paste0("pct:", EDGE_PCTS)))
  st <- lapply(names(thr), function(id) {
    A <- if (startsWith(id, "rho:")) C > thr[[id]] else C >= thr[[id]]
    deg <- rowSums(A); wdeg <- rowSums(C * A)
    list(thr = thr[[id]], deg = deg, ord = order(-deg, -wdeg, runif(n)), n_edges = sum(A) / 2, n_conn = sum(deg >= 1))
  })
  names(st) <- names(thr)
  list(st = st, genes = rownames(C), n = n, npairs = npairs)
}
# number of genes a degree rule::selects, and which genes
K_for <- function(ns, edge_id, deg_type, deg_value) {
  s <- ns$st[[edge_id]]
  if (deg_type == "min") sum(s$deg >= deg_value) else min(ceiling(deg_value / 100 * ns$n), s$n_conn)
}
sel_for <- function(ns, edge_id, K) ns$genes[ns$st[[edge_id]]$ord[seq_len(K)]]

# --------------------------------------------------------------------------
# one network variant
# --------------------------------------------------------------------------
run_variant <- function(variant) {
  cat(sprintf("\n################################################################\n# VARIANT: %s\n################################################################\n", variant_desc[[variant]]))
  net_cells <- if (variant == "induced") score_cells else colnames(counts)
  cn <- counts[, net_cells, drop = FALSE]
  cand_all <- unique(intersect(targetscan_all$ensembl_gene_id, rownames(cn)))
  cv <- apply(cn[cand_all, , drop = FALSE], 1, var)
  cand_ids <- cand_all[is.finite(cv) & cv > 0]; n_cand <- length(cand_ids)
  cat(sprintf("Network cells: %d | candidates: %d | candidate pairs: %s\n", length(net_cells), n_cand, format(n_cand * (n_cand - 1) / 2, big.mark = ",")))

  C_obs <- spearman_full_cpp(cn[cand_ids, , drop = FALSE]); rownames(C_obs) <- colnames(C_obs) <- cand_ids; diag(C_obs) <- 0
  obs_ns <- network_stats(C_obs)

  # ---- how specific numbers and percentages map onto each other in THIS network ----
  mapping <- bind_rows(lapply(names(obs_ns$st), function(id) {
    s <- obs_ns$st[[id]]; tibble(rule = ifelse(startsWith(id, "rho:"), sprintf("rho > %s", sub("rho:", "", id)), sprintf("top %s%% of pairs", sub("pct:", "", id))),
                                 rho_cutoff = s$thr, n_edges = s$n_edges, pct_of_pairs = 100 * s$n_edges / obs_ns$npairs, genes_with_an_edge = s$n_conn)
  }))
  cat("\nEdge rules in this network (specific numbers <-> percentages):\n"); print(as.data.frame(mapping %>% mutate(rho_cutoff = round(rho_cutoff, 3), pct_of_pairs = round(pct_of_pairs, 4))), row.names = FALSE)

  # ---- yardsticks that depend only on K (cached) ----
  ts_ranked <- names(sort(ts_score_gene[cand_ids]))              # strongest TargetScan candidates first
  ref_cache <- new.env()
  ref_for_K <- function(K) {
    key <- as.character(K); if (!is.null(ref_cache[[key]])) return(ref_cache[[key]])
    nl <- replicate(N_NULL_GRID, rho_fast(sample(cand_ids, K)))
    out <- c(ts_topK = rho_fast(ts_ranked[seq_len(K)]), null_median = median(nl, na.rm = TRUE), null_p95 = unname(quantile(nl, 0.95, na.rm = TRUE)))
    assign(key, out, envir = ref_cache); out
  }

  # ---- observed candidates over the whole grid ----
  # Identical gene sets reached by different cutoffs get the same set_id. Sets are stored in a list and compared directly
  # (an earlier version used the pasted gene list as an environment variable NAME, which fails above 10,000 bytes, i.e. ~600 genes).
  set_env <- new.env(); set_env$store <- list()
  find_set_id <- function(sel) {
    ss <- sort(sel)
    for (j in seq_along(set_env$store)) if (length(set_env$store[[j]]) == length(ss) && identical(set_env$store[[j]], ss)) return(j)
    set_env$store[[length(set_env$store) + 1L]] <- ss
    length(set_env$store)
  }
  rows <- lapply(seq_len(nrow(gdf)), function(i) {
    g <- gdf[i, ]; s <- obs_ns$st[[g$edge_id]]
    K <- K_for(obs_ns, g$edge_id, g$deg_type, g$deg_value)
    sel <- if (K >= MIN_K) sel_for(obs_ns, g$edge_id, K) else character(0)
    acc <- if (K >= MIN_K) rho_fast(sel) else NA_real_
    set_id <- if (K >= MIN_K) find_set_id(sel) else NA_integer_
    ref <- if (K >= MIN_K) ref_for_K(K) else c(ts_topK = NA_real_, null_median = NA_real_, null_p95 = NA_real_)
    tibble(grid_id = g$grid_id, set_id = set_id, K = K, K_pct = 100 * K / n_cand, n_edges = s$n_edges, edges_pct = 100 * s$n_edges / obs_ns$npairs,
           rho_cutoff = s$thr, min_deg_selected = if (K > 0) as.numeric(s$deg[s$ord[K]]) else NA_real_,
           accuracy = acc, ts_topK = ref[["ts_topK"]], null_median = ref[["null_median"]], null_p95 = ref[["null_p95"]])
  })
  tab <- left_join(gdf, bind_rows(rows), by = "grid_id") %>% mutate(variant = variant)

  # ---- the same rule on expression-matched NON-target pools ----
  var_n <- apply(cn, 1, var); mean_n <- rowMeans(cn)
  bin_n <- cut(rank(mean_n, ties.method = "average") / length(mean_n), breaks = seq(0, 1, length.out = N_EXPR_BINS + 1), include.lowest = TRUE, labels = FALSE)
  names(bin_n) <- rownames(cn)
  nontarget_ids <- setdiff(rownames(cn), targetscan_all$ensembl_gene_id)
  nontarget_ids <- nontarget_ids[is.finite(var_n[nontarget_ids]) & var_n[nontarget_ids] > 0]
  need <- table(factor(bin_n[cand_ids], levels = seq_len(N_EXPR_BINS)))
  pool_by_bin <- split(nontarget_ids, factor(bin_n[nontarget_ids], levels = seq_len(N_EXPR_BINS)))
  draw_matched_pool <- function() unlist(lapply(seq_len(N_EXPR_BINS), function(b) {
    n_need <- as.integer(need[b]); avail <- pool_by_bin[[b]]
    if (n_need == 0 || length(avail) == 0) return(NULL)
    sample(avail, min(n_need, length(avail)))
  }), use.names = FALSE)

  cat(sprintf("\nRunning %d expression-matched random pools through all %d settings (same K as the candidates at each setting)...\n", N_POOLS_GRID, nrow(gdf)))
  pool_rho <- matrix(NA_real_, nrow(gdf), N_POOLS_GRID)
  for (b in seq_len(N_POOLS_GRID)) {
    pool <- draw_matched_pool()
    Cp <- spearman_full_cpp(cn[pool, , drop = FALSE]); rownames(Cp) <- colnames(Cp) <- pool; diag(Cp) <- 0
    ns <- network_stats(Cp)
    pool_rho[, b] <- vapply(seq_len(nrow(gdf)), function(i) if (tab$K[i] >= MIN_K) rho_fast(sel_for(ns, gdf$edge_id[i], tab$K[i])) else NA_real_, numeric(1))
    if (b %% 10 == 0) cat(sprintf("  ... %d/%d pools\n", b, N_POOLS_GRID))
  }
  n_valid <- rowSums(!is.na(pool_rho))
  tab$pool_median <- apply(pool_rho, 1, median, na.rm = TRUE); tab$pool_median[n_valid == 0] <- NA
  tab$pool_p95 <- apply(pool_rho, 1, function(x) if (all(is.na(x))) NA_real_ else unname(quantile(x, 0.95, na.rm = TRUE)))
  tab$p_emp <- ifelse(is.na(tab$accuracy) | n_valid == 0, NA_real_, (1 + rowSums(pool_rho >= tab$accuracy, na.rm = TRUE)) / (1 + n_valid))
  tab <- tab %>% group_by(set_id) %>% mutate(n_settings_same_set = ifelse(is.na(set_id), NA_integer_, n())) %>% ungroup() %>%
    mutate(advantage_vs_ts = accuracy - ts_topK, margin_vs_pool_p95 = accuracy - pool_p95)

  # ---- console views ----
  fmt_matrix <- function(df, et, dt, fun, title) {
    d <- df %>%dplyr::filter(edge_type == et, deg_type == dt) %>% mutate(cell = fun(.))
    m <- d %>%dplyr::select(edge_value, deg_value, cell) %>% pivot_wider(names_from = deg_value, values_from = cell) %>% arrange(edge_value)
    rn <- if (et == "rho") sprintf("rho>%.2f", m$edge_value) else sprintf("top%g%%pairs", m$edge_value)
    out <- as.data.frame(m[, -1, drop = FALSE]); rownames(out) <- rn
    colnames(out) <- if (dt == "min") paste0(">=", colnames(out), " edges") else paste0("top", colnames(out), "%")
    cat("\n", title, "\n", sep = ""); print(out, right = TRUE)
  }
  f_acc <- function(d) ifelse(is.na(d$accuracy), "-", sprintf("%.2f (%d)", d$accuracy, d$K))
  f_adv <- function(d) ifelse(is.na(d$margin_vs_pool_p95), "-", sprintf("%+.2f%s", d$margin_vs_pool_p95, ifelse(!is.na(d$p_emp) & d$p_emp <= 0.05, "*", "")))
  cat(sprintf("\n=== %s ===\nCells show: accuracy (K = genes::selected). 'margin' = accuracy minus the 95th percentile of the SAME rule on expression-matched random pools; * = observed beats >= 95%% of pools (p <= 0.05).\n", variant_desc[[variant]]))
  fmt_matrix(tab, "rho", "min", f_acc, "ABSOLUTE rho cutoff x ABSOLUTE minimum edges per gene -- accuracy (K):")
  fmt_matrix(tab, "rho", "min", f_adv, "ABSOLUTE rho cutoff x ABSOLUTE minimum edges per gene -- margin over random pools:")
  fmt_matrix(tab, "pct", "pct", f_acc, "PERCENTAGE edges (top % of pairs) x PERCENTAGE genes (top % by connectivity) -- accuracy (K):")
  fmt_matrix(tab, "pct", "pct", f_adv, "PERCENTAGE edges x PERCENTAGE genes -- margin over random pools:")

  ok <- tab %>%dplyr::filter(!is.na(accuracy))
  preds <- c(`rho cutoff` = "rho_cutoff", `min edges ofdplyr::selected genes` = "min_deg_selected", `number of edges` = "n_edges", `edges (% of pairs)` = "edges_pct", `genes selected K` = "K", `genes selected (%)` = "K_pct")
  outs  <- c(accuracy = "accuracy", `margin vs TargetScan top-K` = "advantage_vs_ts", `margin vs pools' p95` = "margin_vs_pool_p95")
  cors <- sapply(outs, function(o) sapply(preds, function(p) suppressWarnings(cor(ok[[p]], ok[[o]], method = "spearman", use = "complete.obs"))))
  rownames(cors) <- names(preds)
  cat(sprintf("\nSpearman correlation across the %d grid settings with K >= %d, whichselect only %d distinct gene sets (DESCRIPTIVE: settings are nested and partly redundant, so no p-values):\n", nrow(ok), MIN_K, length(unique(ok$set_id))))
  print(round(cors, 2))

  best <- ok %>%dplyr::filter(K >= 5, !is.na(margin_vs_pool_p95)) %>% arrange(desc(margin_vs_pool_p95)) %>% distinct(set_id, .keep_all = TRUE) %>% head(8) %>%
    transmute(edge_rule = ifelse(edge_type == "rho", sprintf("rho > %g", edge_value), sprintf("top %g%% pairs", edge_value)),
              degree_rule = ifelse(deg_type == "min", sprintf(">= %g edges", deg_value), sprintf("top %g%% genes", deg_value)),
              also_reached_by = n_settings_same_set - 1L, K, edges = n_edges, accuracy = round(accuracy, 3), TargetScan_topK = round(ts_topK, 3),
              pool_p95 = round(pool_p95, 3), margin = round(margin_vs_pool_p95, 3), p = round(p_emp, 3))
  cat("\nDistinct gene sets with the largest margin over random pools (K >= 5; 'also_reached_by' = other settings giving the identical set).\nNOTE: picked on the same cells they are scored on -> optimistic:\n"); print(as.data.frame(best), row.names = FALSE)

  list(variant = variant, tab = tab, mapping = mapping, cors = cors, pool_rho = pool_rho, n_cand = n_cand, npairs = obs_ns$npairs, n_net_cells = length(net_cells))
}

results <- lapply(VARIANTS, run_variant); names(results) <- VARIANTS
tab_all <- bind_rows(lapply(results, `[[`, "tab"))

# anchors, for reference
cat(sprintf("\nAnchors (induced cells vs fluorescence): all %d candidates %.3f | TargetScan top 200 (%d genes) %.3f\n", length(cand_global), rho_fast(cand_global), length(top200_genes), rho_fast(top200_genes)))

# --------------------------------------------------------------------------
# plots
# --------------------------------------------------------------------------
axis_lab <- function(et, dt) list(
  y = if (et == "rho") "rho cutoff (edge if rho > value)" else "edges = top x% of all candidate pairs",
  x = if (dt == "min") "minimum number of edges per gene (>=)" else "genes kept = top x% by connectivity")
heat <- function(df, et, dt, kind) {
  d <- df %>%dplyr::filter(edge_type == et, deg_type == dt) %>% mutate(x = factor(deg_value), y = factor(edge_value), variant_lab = variant_desc[variant])
  al <- axis_lab(et, dt)
  if (kind == "acc") {
    d <- d %>% mutate(lab = ifelse(is.na(accuracy), "-", sprintf("%.2f\nK=%d", accuracy, K)))
    p <- ggplot(d, aes(x, y, fill = accuracy)) + geom_tile(colour = "white") + geom_text(aes(label = lab), size = 2.6, lineheight = 0.9) +
      scale_fill_viridis_c(option = "C", na.value = "grey90", name = "rho vs\nfluorescence") + labs(title = "Accuracy against GFP fluorescence (induced cells)", subtitle = "tile text: accuracy and K = number of genesdplyr::selected")
  } else {
    lim <- max(abs(d$margin_vs_pool_p95), na.rm = TRUE)
    d <- d %>% mutate(lab = ifelse(is.na(margin_vs_pool_p95), "-", sprintf("%+.2f%s", margin_vs_pool_p95, ifelse(!is.na(p_emp) & p_emp <= 0.05, "*", ""))))
    p <- ggplot(d, aes(x, y, fill = margin_vs_pool_p95)) + geom_tile(colour = "white") + geom_text(aes(label = lab), size = 2.7) +
      scale_fill_gradient2(low = "#b2182b", mid = "white", high = "#2166ac", midpoint = 0, limits = c(-lim, lim), na.value = "grey90", name = "accuracy minus\npools' 95th pct") +
      labs(title = "Margin over the same rule applied to expression-matched random pools", subtitle = "margin = accuracy minus the 95th percentile of the pools; * = beats >= 95% of pools; blue = specific to the candidates, red = within generic structure")
  }
  p + facet_wrap(~variant_lab) + labs(x = al$x, y = al$y) + theme_bw(base_size = 11) + theme(panel.grid = element_blank(), plot.subtitle = element_text(size = 9))
}
for (nm in c("absolute", "percentage")) {
  et <- if (nm == "absolute") "rho" else "pct"; dt <- if (nm == "absolute") "min" else "pct"
  pdf(file.path(OUT_DIR, sprintf("cutoff_benchmark_%s.pdf", nm)), width = 14, height = 6.5); print(heat(tab_all, et, dt, "acc")); print(heat(tab_all, et, dt, "adv")); invisible(dev.off())
}

size_lab <- c(K = "genesdplyr::selected (K)", K_pct = "genesdplyr::selected (% of candidates)", n_edges = "edges (number)", edges_pct = "edges (% of all pairs)")
long <- tab_all %>%dplyr::filter(!is.na(accuracy)) %>%dplyr::select(variant, accuracy, ts_topK, K, K_pct, n_edges, edges_pct) %>%
  pivot_longer(c(K, K_pct, n_edges, edges_pct), names_to = "x_var", values_to = "x") %>%dplyr::filter(x > 0) %>% mutate(x_var = factor(size_lab[x_var], levels = size_lab))
# Spearman correlations go into the facet titles (annotations inside the panels collided with the data)
ann <- long %>% group_by(variant, x_var) %>% summarise(r = suppressWarnings(cor(x, accuracy, method = "spearman")), .groups = "drop")
facet_lab <- ann %>% group_by(x_var) %>% summarise(x_lab = sprintf("%s\nSpearman with accuracy: %s", as.character(x_var[1]), paste(sprintf("%s %.2f", variant, r), collapse = " | ")), .groups = "drop") %>%
  arrange(match(as.character(x_var), size_lab))
long <- long %>% left_join(facet_lab, by = "x_var") %>% mutate(x_lab = factor(x_lab, levels = facet_lab$x_lab))
p_size <- ggplot(long, aes(x, accuracy)) +
  geom_point(aes(y = ts_topK), colour = "grey65", shape = 1, size = 1.3) +
  geom_point(aes(colour = variant), alpha = 0.7, size = 1.5) +
  facet_wrap(~x_lab, scales = "free_x") + scale_x_log10() + theme_bw(base_size = 12) + theme(legend.position = "bottom") +
  labs(x = NULL, y = "rho vs GFP fluorescence (induced cells)", colour = "network built on",
       title = "Accuracy versus how many genes / edges a cutoff keeps",
       subtitle = "coloured = co-expression selection, every grid setting; grey circles = TargetScan's own top-K at the same K (same x)")
ggsave(file.path(OUT_DIR, "cutoff_benchmark_vs_size.pdf"), p_size, width = 11, height = 8)

write.csv(tab_all %>%dplyr::select(-edge_id, -deg_id), file.path(OUT_DIR, "cutoff_benchmark_table.csv"), row.names = FALSE)
saveRDS(list(variants = results, table = tab_all, grid = list(rho = RHO_CUTOFFS, min_degree = MIN_DEGREES, edge_pct = EDGE_PCTS, degree_pct = DEG_PCTS),
             params = list(n_pools = N_POOLS_GRID, n_null = N_NULL_GRID, min_K = MIN_K, n_scoring_cells = length(score_cells), used_cpp = USE_CPP)),
        file.path(OUT_DIR, "res_coexpression_cutoff_benchmark.rds"))

cat("\n================ How to read this ================\n",
    "* Specific numbers vs percentages: the mapping table shows what each percentage means as a rho cutoff and edge count in each network. Percentage rules\n",
    "  keep the same number of edges / genes whatever the network, so they are the fair way to compare the induced and all-cell networks.\n",
    "* Accuracy alone mixes the filter's own effect with how many genes it keeps. The yardsticks separate them: TargetScan top-K (same K), random K-subsets of\n",
    "  the candidates, and the SAME rule on expression-matched random pools (same K). A positive 'margin' (accuracy above the pools' 95th percentile, marked *) means the setting finds something the\n",
    "  candidates have and generic co-expression does not; it does NOT show that this is miR-124 activity (3'UTR length, tissue modules remain possible).\n",
    "* The correlations across grid settings are descriptive. Strong negative correlations with K just say that stricter cutoffs keep fewer genes.\n",
    "* Do not quote the best setting's accuracy as the method's accuracy: it was chosen on these same cells.\n", sep = "")
cat(sprintf("\nSaved cutoff_benchmark_table.csv, three PDFs and res_coexpression_cutoff_benchmark.rds to %s\n", OUT_DIR))
