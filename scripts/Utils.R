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


