library(tidyverse)


calculate_activity <- function(counts, targets) {
  targets <- intersect(targets, rownames(counts))
  cell_sums <- colSums(counts[targets,,drop=FALSE])
  cell_sums_norm <- cell_sums / mean(cell_sums)
  -log2(cell_sums_norm) + 1
}


get_top_miRNAs <- function(miRNA_mat, N = 5) {
  miRNA_mat |>
    rowSums() |>
    sort(decreasing = TRUE) |>
    head(N) |>
    names()
}

mirna_in_db <- function(mir, db_key) {
  mir %in% keys(miRNAtap.db, keytype = db_key)
}

get_valid_databases_for_mir <- function(mir, db_list) {
  present <- vapply(
    db_list,
    function(db) {
      keytype <- db_key_map[[db]]
      mir %in% keys(miRNAtap.db, keytype = keytype)
    },
    logical(1)
  )
  db_list[present]
}

db_key_map <- list(
  diana      = "DIANA_MAPPED_HSA",
  miranda    = "MIRANDA_MAPPED_HSA",
  mirdb      = "MIRDB_MAPPED_HSA",
  pictar     = "PICTAR_MAPPED_HSA",
  targetscan = "TARGETSCAN_MAPPED_HSA",
  clash      = "CLASH_HEWLAK_MAPPED_HSA"
)

# mir <- "hsa-miR-124-3p"
# target_databases <- c('targetscan','diana','mirdb','pictar','miranda')

get_targets_miRNAtap <- function(
    mir,
    target_databases = c('targetscan', 'diana', 'mirdb'),
    species = 'hsa',
    ranking_method = 'geom',
    min_databases = 2,
    min_targets = 20   # minimum acceptable target count
) {
  message("→ Getting targets for: ", mir)
  
  #---------------------------------------------------------
  # STEP 1 — Predict targets (safe)
  #---------------------------------------------------------
  preds <- tryCatch({
    getPredictedTargets(
      mir, species = species,
      method = ranking_method,
      sources = target_databases,
      min_src = min_databases
    ) |> as.data.frame()
  }, error = function(e) {
    warning("⚠️ miRNAtap failed for ", mir, ": ", e$message)
    return(data.frame())  # return an empty data.frame
  })
  
  # No predictions
  if (nrow(preds) == 0) {
    warning("⚠️ No predicted targets found for ", mir)
    return(tibble(
      entrez_id = character(),
      external_gene_name = character(),
      ensembl_gene_id = character(),
      source_miRNA = mir
    ))
  }
  
  #---------------------------------------------------------
  # STEP 2 — Annotate genes via biomaRt (safe)
  #---------------------------------------------------------
  genes <- tryCatch({
    getBM(
      attributes = c('ensembl_gene_id', 'entrezgene_id', 'external_gene_name'),
      filters = "entrezgene_id",
      values = rownames(preds),
      mart = ensembl_hs
    ) |> mutate(entrezgene_id = as.character(entrezgene_id))
  }, error = function(e) {
    warning("⚠️ biomaRt lookup failed for ", mir, ": ", e$message)
    return(data.frame())
  })
  
  if (nrow(genes) == 0) {
    warning("⚠️ biomaRt returned no annotations for ", mir)
    return(tibble(
      entrez_id = rownames(preds),
      external_gene_name = NA,
      ensembl_gene_id = NA,
      source_miRNA = mir
    ))
  }
  
  #---------------------------------------------------------
  # STEP 3 — Merge predictions with annotations
  #---------------------------------------------------------
  targets <- preds |>
    rownames_to_column("entrez_id") |>
    inner_join(genes, join_by(entrez_id == entrezgene_id)) |>
    mutate(source_miRNA = mir)
  
  #---------------------------------------------------------
  # STEP 4 — Check for "too few" targets
  #---------------------------------------------------------
  if (nrow(targets) < min_targets) {
    warning("⚠️ Only ", nrow(targets), " targets for ", mir,
            " (< ", min_targets, "). Proceed carefully.")
  }
  
  return(targets)
}


make_pairs_by_gene_name <- function(targets, N) {
  top_targets <- targets |>
    slice_head(n = N)
  pairs_df <- as.data.frame(t(combn(top_targets$external_gene_name , 2)))
  colnames(pairs_df) <- c("gene1", "gene2")
  return(pairs_df)
}

gene_filter <- function(x) {
  nonzero <- sum(x > 0, na.rm = TRUE)
  zero_prop <- mean(x == 0, na.rm = TRUE)
  nonmiss <- sum(!is.na(x))
  v <- var(x, na.rm = TRUE)
  
  (nonzero >= min_nonzero) &&
    (zero_prop <= max_zero_prop) &&
    (nonmiss >= min_nonmissing) &&
    (!is.na(v) && v > min_var)
}



calculate_pairwise_cors <- function(pairs_df, mRNA){
  tic()
  rho <- furrr::future_map2_dbl(
    pairs_df$gene1, pairs_df$gene2,
    ~ safe_spearman(mRNA[.x, ], mRNA[.y, ])
  ) 
  toc()
  pairwise_rho <- pairs_df |>
    cbind(rho)
}



filter_network <- function(cor_results, rho_cutoff,degree_N ) {
  cor_results_filtered <- cor_results |>
    filter(rho > rho_cutoff & rho < 0.99) |> 
    distinct_all()
  g <- igraph::graph_from_data_frame(cor_results_filtered, directed = FALSE)
  # Degrees
  deg <- igraph::degree(g, mode = "all")
  nodes_to_keep <- names(deg[deg >= degree_N])
}


safe_spearman <- function(x, y) {
  
  x <- as.numeric(x)
  y <- as.numeric(y)
  
  # Remove cases with NA or Inf
  df <- tibble(x = x, y = y) %>% 
    filter(is.finite(x), is.finite(y))
  
  # Need at least 3 observations
  if (nrow(df) < 3) return(NA_real_)
  
  # Need some variation
  if (sd(df$x) == 0 || sd(df$y) == 0) return(NA_real_)
  
  # Try cor.test safely
  out <- suppressWarnings(try(cor.test(df$x, df$y, method = "spearman"), silent = TRUE))
  
  if (inherits(out, "try-error")) return(NA_real_)
  
  out$estimate
}

tic <- function() {
  tic_start <<- base::Sys.time()
}

toc <- function() {
  dt <- base::difftime(base::Sys.time(), tic_start)
  dt <- round(dt, digits = 1L)
  message(paste(format(dt), "since tic()"))
}


########################

# run_mhg_for_column <- function(column_name, counts_df, gene_set, n_max = 1000) {
#   sorted_genes <- counts_df |>
#     rownames_to_column("genes") |>
#     select(genes, all_of(column_name)) |>
#     arrange(.data[[column_name]]) |>
#     pull(genes)
#   
#   v <- as.integer(sorted_genes %in% gene_set)
#   
#   if (sum(v) == 0 || all(is.na(v))) {
#     return(data.frame(
#       column = column_name,
#       stat = NA,
#       cutoff = NA,
#       pval = NA
#     ))
#   }
#   
#   res <- mHG.test(v, n_max = n_max)
#   
#   data.frame(
#     column = column_name,
#     stat = res$statistic,
#     pval = res$p.value
#   )
# }


calculate_activity_mHG <- function(counts, targets){
  gene_set <- targets$ensembl_gene_id
  # Run for all columns (except rownames)
  mhg_results <- future_map(
    colnames(counts),
    run_mhg_for_column,
    counts_df = counts,
    gene_set = gene_set,
    n_max = 1000
  )
  # Combine into a single data frame
  mhg_df <- do.call(rbind, mhg_results)
}




## ============================================================================
## GFP ground truth: mRNA + FACS fluorescence, all cells vs Dox-induced cells
## ----------------------------------------------------------------------------
## Shared by every HEK_SS3 comparison script so all tools are scored on the
## SAME cells with the SAME rules.
##
## QC rule (unchanged from the original scripts): drop cells at or above the
## 99th percentile of scran-normalized GFP mRNA, computed across all cells.
## That same QC-passing set is used for BOTH readouts, so any difference
## between mRNA and fluorescence rho reflects the readout, not a different
## set of cells.
##
## "induced" = Dox_Group == "Dox Induced" (0.01, 0.1, 1 ug/ml Dox).
## "No Dox" covers the 0-Dox wells and the Control wells (P23/P24).
## ============================================================================

GFP_READOUTS <- c(mRNA = "gfp_mrna", fluorescence = "gfp_facs")
CELL_SETS    <- c("all", "induced")

load_gfp_truth <- function(processed_dir, upper_pct = 0.99) {
  path <- file.path(processed_dir, "gfp_correlation_data.csv")
  if (!file.exists(path)) stop("GFP ground-truth table not found: ", path)
  md <- read.csv(path, stringsAsFactors = FALSE)

  needed  <- c("cell_id", "GFP_normalized", "GFP_fluorescence", "Dox_Group", "Dox_Concentration")
  missing <- setdiff(needed, colnames(md))
  if (length(missing) > 0) stop("gfp_correlation_data.csv is missing columns: ", paste(missing, collapse = ", "))
  if (anyDuplicated(md$cell_id)) stop("Duplicate cell_id values in gfp_correlation_data.csv")
  unexpected <- setdiff(unique(md$Dox_Group), c("Dox Induced", "No Dox"))
  if (length(unexpected) > 0) stop("Unexpected Dox_Group labels: ", paste(unexpected, collapse = ", "))

  truth <- data.frame(
    Cell_ID  = md$cell_id,
    gfp_mrna = md$GFP_normalized,
    gfp_facs = md$GFP_fluorescence,
    dox      = md$Dox_Concentration,
    induced  = md$Dox_Group == "Dox Induced",
    stringsAsFactors = FALSE
  )
  truth <- truth[is.finite(truth$gfp_mrna), ]

  cutoff <- unname(quantile(truth$gfp_mrna, upper_pct))
  truth$qc_pass <- truth$gfp_mrna < cutoff & is.finite(truth$gfp_facs)
  attr(truth, "gfp_mrna_cutoff") <- cutoff

  cat(sprintf(paste0("GFP truth: %d cells; QC (GFP mRNA < %.0fth pct = %.3f) keeps %d ",
                     "(%d induced, %d uninduced); %d cells have negative FACS values\n"),
              nrow(truth), 100 * upper_pct, cutoff, sum(truth$qc_pass),
              sum(truth$qc_pass & truth$induced), sum(truth$qc_pass & !truth$induced),
              sum(truth$qc_pass & truth$gfp_facs < 0)))
  truth
}

# Run `expr` with a fixed seed without disturbing the caller's RNG stream.
.with_seed <- function(seed, expr) {
  had_seed <- exists(".Random.seed", envir = globalenv(), inherits = FALSE)
  if (had_seed) old_seed <- get(".Random.seed", envir = globalenv())
  on.exit(if (had_seed) assign(".Random.seed", old_seed, envir = globalenv())
          else rm(".Random.seed", envir = globalenv()))
  set.seed(seed)
  expr
}

# Percentile bootstrap CI for Spearman rho (cells resampled with replacement).
# cor(method = "spearman") uses average ranks for ties, like cor.test().
spearman_boot_ci <- function(x, y, n_boot = 1000, seed = 1312, level = 0.95) {
  if (n_boot <= 0 || length(x) < 3) return(c(NA_real_, NA_real_))
  n <- length(x)
  boots <- .with_seed(seed, replicate(n_boot, {
    i <- sample.int(n, n, replace = TRUE)
    suppressWarnings(cor(x[i], y[i], method = "spearman"))
  }))
  a <- (1 - level) / 2
  unname(quantile(boots, c(a, 1 - a), na.rm = TRUE))
}

subset_cells <- function(d, cell_set) {
  switch(cell_set,
         all     = d,
         induced = d[d$induced, , drop = FALSE],
         stop("Unknown cell_set: ", cell_set))
}

# Score one method: 2 readouts x 2 cell sets -> 4 rows.
# `scores` must have a Cell_ID column and the score column `score_col`
# (sign convention: higher = more miRNA activity).
evaluate_score_vs_gfp <- function(scores, score_col, method_label, truth,
                                  n_boot = 1000, seed = 1312) {
  if (!all(c("Cell_ID", score_col) %in% colnames(scores)))
    stop(sprintf("[%s] scores must contain columns 'Cell_ID' and '%s'", method_label, score_col))

  s <- data.frame(Cell_ID = scores$Cell_ID, score = scores[[score_col]], stringsAsFactors = FALSE)
  d <- merge(s, truth[truth$qc_pass, ], by = "Cell_ID")
  d <- d[is.finite(d$score), ]
  if (nrow(d) == 0) stop(sprintf("[%s] no cells left after matching to GFP truth -- check Cell_ID format", method_label))

  out <- list()
  for (cs in CELL_SETS) {
    sub <- subset_cells(d, cs)
    for (ro in names(GFP_READOUTS)) {
      y  <- sub[[GFP_READOUTS[[ro]]]]
      ct <- suppressWarnings(cor.test(sub$score, y, method = "spearman", exact = FALSE))
      ci <- spearman_boot_ci(sub$score, y, n_boot = n_boot, seed = seed)
      out[[length(out) + 1]] <- tibble::tibble(
        method = method_label, readout = ro, cell_set = cs,
        rho = unname(ct$estimate), ci_low = ci[1], ci_high = ci[2],
        p_value = ct$p.value, n = nrow(sub)
      )
    }
  }
  res <- dplyr::bind_rows(out)
  for (k in seq_len(nrow(res))) {
    cat(sprintf("  %-38s %-12s %-8s rho = %.3f [%.3f, %.3f]  n = %d\n",
                res$method[k], res$readout[k], res$cell_set[k],
                res$rho[k], res$ci_low[k], res$ci_high[k], res$n[k]))
  }
  res
}

# Paired comparison of two methods on the SAME cells: delta = rho_a - rho_b,
# with a bootstrap CI from resampling cells jointly for both methods.
# Also reports rho between the two scores themselves (needed to judge
# whether two methods that both track GFP are really capturing the same thing).
paired_rho_difference <- function(scores_a, col_a, label_a,
                                  scores_b, col_b, label_b,
                                  truth, n_boot = 1000, seed = 1312) {
  a <- data.frame(Cell_ID = scores_a$Cell_ID, score_a = scores_a[[col_a]], stringsAsFactors = FALSE)
  b <- data.frame(Cell_ID = scores_b$Cell_ID, score_b = scores_b[[col_b]], stringsAsFactors = FALSE)
  d <- merge(merge(a, b, by = "Cell_ID"), truth[truth$qc_pass, ], by = "Cell_ID")
  d <- d[is.finite(d$score_a) & is.finite(d$score_b), ]

  out <- list()
  for (cs in CELL_SETS) {
    sub <- subset_cells(d, cs)
    n <- nrow(sub)
    for (ro in names(GFP_READOUTS)) {
      y <- sub[[GFP_READOUTS[[ro]]]]
      rho_a <- cor(sub$score_a, y, method = "spearman")
      rho_b <- cor(sub$score_b, y, method = "spearman")
      deltas <- .with_seed(seed, replicate(n_boot, {
        i <- sample.int(n, n, replace = TRUE)
        cor(sub$score_a[i], y[i], method = "spearman") - cor(sub$score_b[i], y[i], method = "spearman")
      }))
      ci <- unname(quantile(deltas, c(0.025, 0.975), na.rm = TRUE))
      out[[length(out) + 1]] <- tibble::tibble(
        method_a = label_a, method_b = label_b, readout = ro, cell_set = cs,
        rho_a = rho_a, rho_b = rho_b, delta = rho_a - rho_b,
        delta_ci_low = ci[1], delta_ci_high = ci[2],
        rho_between_methods = cor(sub$score_a, sub$score_b, method = "spearman"),
        n = n
      )
    }
  }
  res <- dplyr::bind_rows(out)
  for (k in seq_len(nrow(res))) {
    cat(sprintf("  %s - %s | %-12s %-8s delta = %+.3f [%+.3f, %+.3f]  rho(a,b) = %.3f  n = %d\n",
                res$method_a[k], res$method_b[k], res$readout[k], res$cell_set[k],
                res$delta[k], res$delta_ci_low[k], res$delta_ci_high[k],
                res$rho_between_methods[k], res$n[k]))
  }
  res
}

# Forest-style plot of stratified results: rows = methods, panels = readout x cell set.
plot_stratified <- function(stratified, title = "miR-124 activity vs GFP (Spearman rho, 95% bootstrap CI)") {
  d <- stratified |>
    dplyr::mutate(
      readout  = factor(readout, levels = names(GFP_READOUTS),
                        labels = c("GFP mRNA (scran)", "GFP fluorescence (FACS)")),
      cell_set = factor(cell_set, levels = CELL_SETS,
                        labels = c("All cells", "Dox-induced cells only")),
      is_ours  = grepl("^Our method", method)
    )
  order_by <- d |> dplyr::filter(readout == "GFP fluorescence (FACS)", cell_set == "Dox-induced cells only")
  lvl <- if (nrow(order_by) > 0) order_by$method[order(order_by$rho)] else unique(d$method)
  d$method <- factor(d$method, levels = unique(c(lvl, unique(d$method))))

  ggplot2::ggplot(d, ggplot2::aes(x = rho, y = method, colour = is_ours)) +
    ggplot2::geom_vline(xintercept = 0, linetype = "dashed", colour = "grey60") +
    ggplot2::geom_errorbarh(ggplot2::aes(xmin = ci_low, xmax = ci_high), height = 0.25) +
    ggplot2::geom_point(size = 2.5) +
    ggplot2::geom_text(ggplot2::aes(label = sprintf("%.2f", rho)), vjust = -0.9, size = 3, show.legend = FALSE) +
    ggplot2::facet_grid(cell_set ~ readout) +
    ggplot2::scale_colour_manual(values = c(`TRUE` = "#b2182b", `FALSE` = "#2166ac"), guide = "none") +
    ggplot2::coord_cartesian(xlim = c(min(0, min(d$ci_low, na.rm = TRUE)), 1)) +
    ggplot2::theme_bw(base_size = 12) +
    ggplot2::theme(panel.grid.minor = ggplot2::element_blank()) +
    ggplot2::labs(x = "Spearman rho with GFP", y = NULL, title = title,
                  subtitle = "Same QC-passing cells in every panel; methods ordered by induced-cell fluorescence rho")
}
