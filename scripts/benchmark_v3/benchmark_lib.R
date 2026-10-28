# ==============================================================================
# Benchmark v3 -- shared functions.
#
# Vocabulary
#   miRNA family  miRNAs sharing a seed (nt 2-8); they share their targets.
#   evidence      for one family, one table per database: which genes it
#                 targets and how strongly (see empty_evidence()).
#   target list   a set of genes built from evidence by one rule ("list type",
#                 e.g. TargetScan_top200 or miRTarBase_high).
#   kind          "benchmark"        list for the benchmarked family
#                 "negative_control" same rule, applied to the control miRNA
#                                    after removing the benchmarked family's
#                                    targets (CONTROL_EXCLUSION)
#                 "random"           random genes, same size, drawn from genes
#                                    that are not targets of the benchmarked
#                                    family (CONTROL_EXCLUSION; N_RANDOM draws)
# ==============================================================================
suppressPackageStartupMessages({ library(tidyverse); library(Matrix) })

# ------------------------------------------------------------------------------
# Small utilities
# ------------------------------------------------------------------------------
read_count_matrix <- function(path, cells_in_rows = FALSE) {
  if (!file.exists(path)) stop("Not found: ", path)
  if (grepl("\\.rds$", path, ignore.case = TRUE)) {
    m <- readRDS(path)
  } else {
    sep <- if (grepl("\\.csv(\\.gz)?$", path, ignore.case = TRUE)) "," else "\t"
    # R write.table() files have one header field fewer than data fields;
    # fread then names the row-name column V1 and warns. That warning is expected.
    df <- withCallingHandlers(data.table::fread(path, sep = sep, data.table = FALSE),
                              warning = function(w) if (grepl("column names but the data has", conditionMessage(w))) invokeRestart("muffleWarning"))
    rownames(df) <- df[[1]]; df[[1]] <- NULL; m <- df
  }
  m <- as.matrix(m)
  if (cells_in_rows) m <- t(m)
  storage.mode(m) <- "double"
  m
}

revcomp_dna <- function(s) {
  s <- toupper(gsub("U", "T", s))
  paste(rev(strsplit(chartr("ACGT", "TGCA", s), "")[[1]]), collapse = "")
}
# The 7-mer target site of a miRNA: reverse complement of its seed (nt 2-8)
seed_of <- function(mature_seq) toupper(gsub("U", "T", substr(mature_seq, 2, 8)))
site_of <- function(mature_seq) revcomp_dna(seed_of(mature_seq))

read_fasta <- function(path) {
  x <- readLines(path)
  is_name <- startsWith(x, ">")
  setNames(x[which(is_name) + 1], sub("^>(\\S+).*$", "\\1", x[is_name]))
}

# Gene identifiers: miReact's 3'UTR table maps transcript, gene ID and symbol
load_utr_table <- function(species) {
  readRDS(file.path(MIREACT_DIR, "seqs", paste0(species, ".utr3.seqs.rds")))
}
symbol_to_gene <- function(symbols, utr, universe) {
  map <- utr |> distinct(gid, gsym) |> filter(gid %in% universe, !is.na(gsym), gsym != "")
  map$gid[match(toupper(symbols), toupper(map$gsym))]               # first match; NA if unknown
}
transcript_to_gene <- function(tx, utr, universe) {
  map <- utr |> distinct(tid, gid) |> mutate(tid = sub("\\.\\d+$", "", tid)) |> filter(gid %in% universe)
  map$gid[match(sub("\\.\\d+$", "", tx), map$tid)]
}

# ------------------------------------------------------------------------------
# Evidence: one table per database for one miRNA family
# ------------------------------------------------------------------------------
empty_evidence <- function() list(
  TargetScan = tibble(gene = character(), score = numeric(),        # cumulative weighted context++ (more negative = stronger)
                      context = numeric(), pct = numeric(),         # total context++ score; aggregate PCT (conservation)
                      n8 = numeric(), n7m8 = numeric(), n7a1 = numeric(),   # sites of each type (conserved + poorly conserved)
                      n_conserved = numeric()),
  TarBase    = tibble(gene = character(), support = numeric()),     # independent records
  miRTarBase = tibble(gene = character(), strong = logical(), very_high = logical()),   # functional evidence only
  mirDIP     = tibble(gene = character(), class = character()),     # best score class over family members
  miRDB      = tibble(gene = character(), score = numeric())        # miRDB target score (higher = stronger)
)

# All genes a database lists as targets of a family. These are removed before
# building that database's negative-control lists and from its random pool.
full_target_list <- function(evidence, database) {
  e <- evidence[[database]]
  if (database == "mirDIP") e <- filter(e, class %in% c("Medium", "High", "Very High"))
  unique(e$gene)
}
# Genes to remove per database for the controls of a benchmarked family
benchmark_exclusion <- function(evidence) {
  per_db <- map(set_names(names(evidence)), \(db) full_target_list(evidence, db))
  if (CONTROL_EXCLUSION == "same_database") return(per_db)
  if (CONTROL_EXCLUSION != "all_databases") stop("CONTROL_EXCLUSION must be 'all_databases' or 'same_database'")
  everything <- unique(unlist(per_db))
  map(per_db, \(x) everything)
}
# Targets of OTHER benchmarked families, removed from a family's own lists so it is
# tested on its own targets (mESC single families; SHARED_TARGET_EXCLUSION)
shared_targets <- function(other_evidence) {
  dbs <- names(empty_evidence())
  per_db <- map(set_names(dbs), \(db) unique(unlist(map(other_evidence, \(e) full_target_list(e, db)))))
  if (SHARED_TARGET_EXCLUSION == "same_database") return(per_db)
  if (SHARED_TARGET_EXCLUSION != "all_databases") stop("SHARED_TARGET_EXCLUSION must be 'same_database' or 'all_databases'")
  everything <- unique(unlist(per_db))
  map(per_db, \(x) everything)
}
remove_targets <- function(evidence, exclude) {
  imap(evidence, function(e, db) filter(e, !gene %in% exclude[[db]]))
}

# Which databases a list type comes from (TarBase_TSsite needs both)
list_type_databases <- function(list_type) {
  if (list_type == "TarBase_TSsite") return(c("TarBase", "TargetScan"))
  sub("_.*$", "", list_type)
}

# ------------------------------------------------------------------------------
# Co-expression filter
#   Spearman correlations among candidate genes; the top COEXP_EDGE_PERCENT of
#   pairs are edges; genes with >= COEXP_MIN_EDGES edges are kept.
#   cell_groups: list of cell-ID vectors. Ranks are z-scored inside each group
#   and pooled, so correlations reflect variation WITHIN groups (one group =
#   plain Spearman over those cells).
# ------------------------------------------------------------------------------
coexpression_filter <- function(genes, expr, cell_groups, detected) {
  candidates <- intersect(genes, detected)
  cells <- unlist(cell_groups)
  M <- expr[candidates, cells, drop = FALSE]
  M <- M[apply(M, 1, var) > 0, , drop = FALSE]
  if (nrow(M) < 10) return(character(0))
  Z <- matrix(0, ncol(M), nrow(M), dimnames = list(colnames(M), rownames(M)))
  for (g in cell_groups) {
    R <- apply(M[, g, drop = FALSE], 1, rank)
    Z[g, ] <- scale(R)
  }
  Z[!is.finite(Z)] <- 0                                   # genes constant within a group
  C <- cor(Z); diag(C) <- NA
  pair_cor <- C[upper.tri(C)]
  n_edges <- max(1, ceiling(COEXP_EDGE_PERCENT / 100 * length(pair_cor)))
  cutoff <- sort(pair_cor, partial = length(pair_cor) - n_edges + 1)[length(pair_cor) - n_edges + 1]
  degree <- rowSums(C >= cutoff, na.rm = TRUE)
  names(degree)[degree >= COEXP_MIN_EDGES]
}

# ------------------------------------------------------------------------------
# From evidence to target lists. The same function builds benchmarked and
# negative-control lists, so both follow identical rules and sizes.
#   coexp(genes) is a function applying coexpression_filter() for the dataset.
# ------------------------------------------------------------------------------
build_target_lists <- function(evidence, coexp, with_mirdip) {
  L <- list()

  ts <- evidence$TargetScan
  ts_ranked <- ts |> arrange(score, gene) |> pull(gene)
  L$TargetScan_all <- ts_ranked
  for (n in TARGETSCAN_TOP_N) L[[paste0("TargetScan_top", n)]] <- head(ts_ranked, n)
  # by best site type (8mer > 7mer-m8 > 7mer-A1), counting conserved and poorly conserved sites
  L$TargetScan_8mer           <- ts$gene[which(ts$n8 >= 1)]
  L$TargetScan_8mer_or_7merm8 <- ts$gene[which(ts$n8 + ts$n7m8 >= 1)]
  L$TargetScan_7merm8_best    <- ts$gene[which(ts$n8 == 0 & ts$n7m8 >= 1)]
  L$TargetScan_7merA1_best    <- ts$gene[which(ts$n8 == 0 & ts$n7m8 == 0 & ts$n7a1 >= 1)]
  # by conservation, and by the alternative scores
  L$TargetScan_conserved      <- ts$gene[which(ts$n_conserved >= 1)]
  by_pct     <- ts |> filter(!is.na(pct)) |> arrange(desc(pct), score, gene) |> pull(gene)
  by_context <- ts |> filter(!is.na(context)) |> arrange(context, gene) |> pull(gene)
  for (n in TARGETSCAN_ALT_TOP_N) {
    L[[paste0("TargetScan_PCT_top", n)]]     <- head(by_pct, n)
    L[[paste0("TargetScan_context_top", n)]] <- head(by_context, n)
  }

  tb <- evidence$TarBase
  L$TarBase_all <- tb$gene
  L[[paste0("TarBase_support", TARBASE_MIN_RECORDS)]] <- tb$gene[tb$support >= TARBASE_MIN_RECORDS]
  L$TarBase_TSsite <- intersect(tb$gene, ts_ranked)        # CLIP-supported AND a predicted seed site

  mt <- evidence$miRTarBase
  L$miRTarBase_all      <- mt$gene
  L$miRTarBase_high     <- mt$gene[mt$strong]
  L$miRTarBase_veryhigh <- mt$gene[mt$very_high]

  if (with_mirdip) {
    dip <- evidence$mirDIP
    L$mirDIP_medium_up <- dip$gene[dip$class %in% c("Medium", "High", "Very High")]
    L$mirDIP_high_up   <- dip$gene[dip$class %in% c("High", "Very High")]
    L$mirDIP_veryhigh  <- dip$gene[dip$class == "Very High"]
  }

  db_ranked <- evidence$miRDB |> arrange(desc(score), gene) |> pull(gene)
  L$miRDB_all <- db_ranked
  for (p in MIRDB_TOP_PERCENT) L[[paste0("miRDB_top", p, "pct")]] <- head(db_ranked, ceiling(p / 100 * length(db_ranked)))

  L$TargetScan_coexp        <- coexp(L$TargetScan_all)
  L$TargetScan_top500_coexp <- coexp(L$TargetScan_top500)
  L$TarBase_coexp           <- coexp(L$TarBase_all)
  L$miRTarBase_coexp        <- coexp(L$miRTarBase_all)
  if (with_mirdip) L$mirDIP_coexp <- coexp(L$mirDIP_medium_up)
  map(L, unique)
}

# Random lists: for each list type, N_RANDOM draws of the benchmarked list's
# size from the universe minus the benchmarked family's targets (exclusion from
# benchmark_exclusion()). Returns NULL for a list type whose pool is smaller
# than the list.
build_random_lists <- function(benchmark_lists, exclusion, universe) {
  imap(benchmark_lists, function(genes, list_type) {
    n <- length(genes)
    excluded <- unique(unlist(exclusion[list_type_databases(list_type)]))
    pool <- setdiff(universe, excluded)
    if (n == 0 || length(pool) < n) return(NULL)
    replicate(N_RANDOM, sample(pool, n), simplify = FALSE)
  })
}

# One family -> rows of the lists table (benchmark, negative control, random)
assemble_list_rows <- function(family, benchmark, control, random) {
  rows <- list(
    tibble(family = family, list_type = names(benchmark), kind = "benchmark", draw = NA_integer_, genes = unname(benchmark)),
    tibble(family = family, list_type = names(control), kind = "negative_control", draw = NA_integer_, genes = unname(control)))
  for (lt in names(random)) if (!is.null(random[[lt]]))
    rows[[length(rows) + 1]] <- tibble(family = family, list_type = lt, kind = "random", draw = seq_along(random[[lt]]), genes = random[[lt]])
  bind_rows(rows)
}

# Random 7-mer motifs (controls for motif-based methods), avoiding given motifs
random_motifs <- function(n, avoid) {
  out <- character(0)
  while (length(out) < n) {
    m <- paste(sample(c("A", "C", "G", "T"), 7, replace = TRUE), collapse = "")
    if (!m %in% c(avoid, out)) out <- c(out, m)
  }
  out
}

# ------------------------------------------------------------------------------
# Saving the standard outputs of a target-list script
# ------------------------------------------------------------------------------
save_target_lists <- function(ds, lists, universe, families, motif_lists, utr, extra = list()) {
  lists <- lists |> mutate(genes = map(genes, \(g) intersect(g, universe)), n_genes = lengths(genes),
                           list_id = sprintf("L%05d", row_number()))
  # genes no benchmarked list of the family contains: background for variant C2
  background <- lists |> filter(kind == "benchmark") |> group_by(family) |>
    summarise(genes = list(setdiff(universe, unlist(genes))), .groups = "drop")
  out <- c(list(dataset = ds$name, expression = ds$expression, universe = universe, families = families,
                lists = lists, motif_lists = motif_lists, background = setNames(background$genes, background$family)), extra)
  saveRDS(out, file.path(ds$out_dir, "target_lists.rds"))

  # symbols for miTEA-HiRes (same gene -> symbol rule as the miTEA count exports)
  g2s <- utr |> distinct(gid, gsym) |> filter(!is.na(gsym), gsym != "") |> filter(!duplicated(gid))
  sym <- lists |> filter(n_genes > 0) |> dplyr::select(list_id, genes) |> unnest(genes) |>
    mutate(gene_symbol = g2s$gsym[match(genes, g2s$gid)]) |> filter(!is.na(gene_symbol)) |>
    distinct(list_id, gene_symbol)
  write_tsv(sym, file.path(ds$out_dir, "target_lists_symbols.tsv"))
  jsonlite::write_json(list(dataset = ds$name, mitea_counts_dir = ds$mitea_counts_dir, out_dir = ds$out_dir),
                       file.path(ds$out_dir, "dataset_info.json"), auto_unbox = TRUE, pretty = TRUE)

  report_list_sizes(lists)
  invisible(out)
}

report_list_sizes <- function(lists) {
  sizes <- lists |> filter(kind != "random") |> dplyr::select(family, list_type, kind, n_genes) |>
    pivot_wider(names_from = kind, values_from = n_genes)
  rnd <- lists |> filter(kind == "random") |> distinct(family, list_type) |> mutate(random = "yes")
  sizes <- sizes |> left_join(rnd, by = c("family", "list_type")) |> mutate(random = replace_na(random, "NO (pool too small)"))
  cat("\n=== Target lists: genes per list (benchmark / negative control) and whether random lists exist ===\n")
  print(as.data.frame(sizes), row.names = FALSE)
  small <- sizes |> filter(benchmark < MIN_LIST_SIZE | negative_control < MIN_LIST_SIZE)
  if (nrow(small)) cat(sprintf("%d list types have a benchmark or control list with < %d genes; they are scored but not ranked.\n",
                               nrow(small), MIN_LIST_SIZE))
}

# ------------------------------------------------------------------------------
# Reading target databases into evidence tables
# ------------------------------------------------------------------------------
# First of the candidate column names present in a file header; clear error otherwise
pick_col <- function(header, candidates, what) {
  hit <- candidates[candidates %in% header]
  if (!length(hit)) stop(sprintf("No column for %s (tried: %s). Columns are: %s", what,
                                 paste(candidates, collapse = ", "), paste(header, collapse = ", ")))
  hit[1]
}
file_header <- function(path) names(data.table::fread(path, nrows = 0))

# TargetScan columns under both naming styles (per-family web downloads / bulk Summary_Counts)
TARGETSCAN_COLUMNS <- list(
  score   = c("Cumulative weighted context++ score"),
  context = c("Total context++ score"),
  pct     = c("Aggregate PCT"),
  c8      = c("Conserved 8mer sites", "Number of conserved 8mer sites"),
  c7m8    = c("Conserved 7mer-m8 sites", "Number of conserved 7mer-m8 sites"),
  c7a1    = c("Conserved 7mer-A1 sites", "Number of conserved 7mer-1a sites", "Number of conserved 7mer-A1 sites"),
  p8      = c("Poorly conserved 8mer sites", "Number of nonconserved 8mer sites"),
  p7m8    = c("Poorly conserved 7mer-m8 sites", "Number of nonconserved 7mer-m8 sites"),
  p7a1    = c("Poorly conserved 7mer-A1 sites", "Number of nonconserved 7mer-1a sites", "Number of nonconserved 7mer-A1 sites"),
  n_conserved = c("Conserved sites total", "Total num conserved sites"),
  # poorly conserved families: one count per site type, no conserved / poorly conserved split
  t8      = c("8mer sites"),
  t7m8    = c("7mer-m8 sites"),
  t7a1    = c("7mer-A1 sites"))
targetscan_column_names <- function(header, source = "TargetScan file") {
  found <- map_chr(TARGETSCAN_COLUMNS, \(cands) { hit <- cands[cands %in% header]; if (length(hit)) hit[1] else NA_character_ })
  if (is.na(found[["score"]])) stop(source, " has no 'Cumulative weighted context++ score' column. Columns: ", paste(header, collapse = ", "))
  split_sites <- !is.na(found[["c8"]]); total_sites <- !is.na(found[["t8"]])
  if (!split_sites && total_sites)
    cat(sprintf("  %s: poorly conserved family format (site counts without a conservation split; no conserved sites)\n", source))
  needed <- c("context", "pct", if (split_sites || !total_sites) c("c8", "c7m8", "c7a1", "p8", "p7m8", "p7a1", "n_conserved") else c("t8", "t7m8", "t7a1"))
  missing <- needed[is.na(found[needed])]
  if (length(missing)) cat(sprintf("*** %s: columns not found (list types using them stay empty): %s\n", source, paste(missing, collapse = ", ")))
  found
}
# raw TargetScan rows -> standard columns; one row per gene (its best-scoring transcript)
standardise_targetscan <- function(raw, found, gene) {
  num <- function(key) if (is.na(found[[key]])) rep(NA_real_, nrow(raw)) else suppressWarnings(as.numeric(raw[[found[[key]]]]))
  # sites of one type = conserved + poorly conserved, or the unsplit count for poorly
  # conserved families; NA only when the file has no site columns at all
  sites <- function(conserved, poorly, total) {
    if (!is.na(found[[conserved]])) return(replace_na(num(conserved), 0) + replace_na(num(poorly), 0))
    if (!is.na(found[[total]])) return(replace_na(num(total), 0))
    rep(NA_real_, nrow(raw))
  }
  tibble(gene = gene, score = num("score"), context = num("context"),
         pct = num("pct"),                                   # "< 0.1" (text) becomes NA: not ranked by PCT
         n8 = sites("c8", "p8", "t8"), n7m8 = sites("c7m8", "p7m8", "t7m8"), n7a1 = sites("c7a1", "p7a1", "t7a1"),
         n_conserved = num("n_conserved")) |>
    filter(!is.na(gene), !is.na(score)) |>
    group_by(gene) |> slice_min(score, n = 1, with_ties = FALSE) |> ungroup()
}

missing_db <- function(path, db) {
  if (file.exists(path)) return(FALSE)
  cat(sprintf("*** %s SKIPPED (file not found: %s) -- its lists will be empty\n", db, path)); TRUE
}
best_class <- function(cls) {                                  # mirDIP class order
  lev <- c("Low", "Medium", "High", "Very High"); lev[suppressWarnings(max(match(cls, lev), na.rm = TRUE))]
}
mirtarbase_evidence <- function(rows) {
  # rows: one per interaction record with columns gene, support, experiments (functional records only)
  rows |> group_by(gene) |>
    summarise(strong = any(support == "Functional MTI"),
              very_high = any(support == "Functional MTI" & grepl("reporter", experiments, ignore.case = TRUE)) &
                          any(support == "Functional MTI" & grepl("western|immunoblot|qRT-PCR|qPCR", experiments, ignore.case = TRUE)),
              .groups = "drop")
}

# Human: TargetScan, TarBase, miRTarBase, mirDIP and miRDB bulk files, read once
# for all families. families: named list, family name -> member miRBase IDs;
# seeds: named vector, family name -> seed (nt 2-8, DNA letters).
load_human_evidence <- function(families, seeds, universe, utr) {
  ev <- setNames(rep(list(empty_evidence()), length(families)), names(families))
  all_members <- unique(unlist(families))

  if (!missing_db(HUMAN_DB$targetscan_summary, "TargetScan")) {
    # In TargetScan 8.0's Summary_Counts the family column holds the SEED (RNA letters)
    h <- file_header(HUMAN_DB$targetscan_summary)
    found <- targetscan_column_names(h, basename(HUMAN_DB$targetscan_summary))
    id_cols <- c(tx = pick_col(h, "Transcript ID", "TargetScan transcript"), symbol = pick_col(h, "Gene Symbol", "TargetScan gene symbol"),
                 seed = pick_col(h, c("miR Family", "miRNA family", "miR family"), "TargetScan family (seed)"),
                 species = pick_col(h, "Species ID", "TargetScan species"))
    ts <- data.table::fread(HUMAN_DB$targetscan_summary, data.table = FALSE, select = unique(c(id_cols, na.omit(found))))
    ts <- ts[ts[[id_cols[["species"]]]] == 9606, , drop = FALSE]
    ts$seed_dna <- toupper(gsub("U", "T", ts[[id_cols[["seed"]]]]))
    ts <- ts[ts$seed_dna %in% seeds, , drop = FALSE]
    gene <- transcript_to_gene(ts[[id_cols[["tx"]]]], utr, universe)
    gene[is.na(gene)] <- symbol_to_gene(ts[[id_cols[["symbol"]]]][is.na(gene)], utr, universe)
    for (f in names(families)) {
      rows <- ts$seed_dna == seeds[[f]]
      ev[[f]]$TargetScan <- standardise_targetscan(ts[rows, , drop = FALSE], found, gene[rows])
    }
  }
  if (!missing_db(HUMAN_DB$tarbase, "TarBase")) {
    h <- file_header(HUMAN_DB$tarbase)
    tb <- data.table::fread(HUMAN_DB$tarbase, data.table = FALSE,
                            select = c(pick_col(h, c("mirna_name", "mirna"), "TarBase miRNA"), pick_col(h, "gene_id", "TarBase gene_id")))
    names(tb) <- c("mirna_name", "gene_id")
    tb <- tb |> filter(mirna_name %in% all_members) |> mutate(gene = sub("\\..*$", "", gene_id)) |> filter(gene %in% universe) |>
      count(mirna_name, gene, name = "records")
    for (f in names(families)) ev[[f]]$TarBase <- tb |> filter(mirna_name %in% families[[f]]) |>
      group_by(gene) |> summarise(support = suppressWarnings(max(records)), .groups = "drop")     # max over members: one CLIP study is not counted per member
  }
  if (!missing_db(HUMAN_DB$mirtarbase, "miRTarBase")) {
    mt <- readr::read_csv(HUMAN_DB$mirtarbase, show_col_types = FALSE)
    h <- names(mt)
    mt <- tibble(mirna = mt[[pick_col(h, "miRNA", "miRTarBase miRNA")]],
                 symbol = mt[[pick_col(h, "Target Gene", "miRTarBase target gene")]],
                 species = mt[[pick_col(h, c("Species (Target Gene)", "Species..Target.Gene."), "miRTarBase target species")]],
                 experiments = mt[[pick_col(h, "Experiments", "miRTarBase experiments")]],
                 support = mt[[pick_col(h, c("Support Type", "Support.Type"), "miRTarBase support type")]]) |>
      filter(species == "Homo sapiens", mirna %in% all_members, support %in% c("Functional MTI", "Functional MTI (Weak)")) |>
      mutate(gene = symbol_to_gene(symbol, utr, universe)) |> filter(!is.na(gene))
    for (f in names(families)) ev[[f]]$miRTarBase <- mirtarbase_evidence(filter(mt, mirna %in% families[[f]]))
  }
  if (!missing_db(HUMAN_DB$mirdip, "mirDIP")) {
    dip <- readr::read_tsv(HUMAN_DB$mirdip, show_col_types = FALSE)
    h <- names(dip)
    dip <- tibble(mirna = dip[[pick_col(h, "MicroRNA", "mirDIP miRNA")]], symbol = dip[[pick_col(h, "Gene Symbol", "mirDIP gene symbol")]],
                  class = dip[[pick_col(h, "Score Class", "mirDIP score class")]]) |>
      filter(mirna %in% all_members) |> mutate(gene = symbol_to_gene(symbol, utr, universe)) |> filter(!is.na(gene))
    not_in_dip <- setdiff(all_members, dip$mirna)
    if (length(not_in_dip)) cat("mirDIP has no entries for:", paste(head(not_in_dip, 20), collapse = ", "),
                                if (length(not_in_dip) > 20) "...", "\n")
    for (f in names(families)) ev[[f]]$mirDIP <- dip |> filter(mirna %in% families[[f]]) |>
      group_by(gene) |> summarise(class = best_class(class), .groups = "drop")
  }
  if (!missing_db(HUMAN_DB$mirdb, "miRDB")) {
    md <- readr::read_tsv(HUMAN_DB$mirdb, col_names = c("mirna", "refseq", "score"), col_types = "ccd") |>
      filter(mirna %in% all_members) |> mutate(refseq = sub("\\.\\d+$", "", refseq))
    rs <- suppressMessages(AnnotationDbi::select(org.Hs.eg.db::org.Hs.eg.db, keys = unique(md$refseq),
                                                 keytype = "REFSEQ", columns = "ENSEMBL"))
    md <- md |> inner_join(rs, by = c("refseq" = "REFSEQ"), relationship = "many-to-many") |>
      filter(ENSEMBL %in% universe) |> rename(gene = ENSEMBL)
    for (f in names(families)) ev[[f]]$miRDB <- md |> filter(mirna %in% families[[f]]) |>
      group_by(gene) |> summarise(score = suppressWarnings(max(score)), .groups = "drop")
  }
  ev
}

# Human miRNA families from TargetScan's miR_Family_Info: every human miRNA
# whose seed (nt 2-8) equals the seed of the given miRNA.
human_family_of <- function(mirna, family_info) {
  seed <- family_info$seed[family_info$mirbase_id == mirna]
  if (!length(seed)) {
    core <- sub("^[a-z]{3}-(miR|let)-([0-9a-z]+).*$", "\\2", mirna)            # e.g. "124" from hsa-miR-124-3p
    similar <- grep(paste0("-", core, "([^0-9]|$)"), family_info$mirbase_id, value = TRUE)
    stop(sprintf("%s is not in miR_Family_Info (human, %d entries read from %s). Similar IDs: %s",
                 mirna, nrow(family_info), HUMAN_DB$targetscan_families,
                 if (length(similar)) paste(head(similar, 10), collapse = ", ") else "none"))
  }
  list(seed = seed, members = sort(family_info$mirbase_id[family_info$seed == seed]))
}
read_family_info <- function() {
  fi <- readr::read_tsv(HUMAN_DB$targetscan_families, show_col_types = FALSE)
  fi <- tibble(mirbase_id = trimws(fi$`MiRBase ID`), species = fi$`Species ID`, ts_family = fi$`miR family`,
               seed = toupper(gsub("U", "T", fi$`Seed+m8`))) |> filter(species == 9606)
  # TargetScan 8.0 gives some miRNAs isoform-specific IDs: hsa-miR-124-3p.1 (canonical
  # seed) and hsa-miR-124-3p.2 (shifted seed). Count tables and the other databases use
  # the plain miRBase name, so keep the canonical entry (.1, or no suffix) under that name.
  isoform <- suppressWarnings(as.integer(sub("^.*\\.(\\d+)$", "\\1", ifelse(grepl("\\.\\d+$", fi$mirbase_id), fi$mirbase_id, NA))))
  fi |> mutate(mirbase_id = sub("\\.\\d+$", "", mirbase_id)) |>
    filter(is.na(isoform) | isoform == 1) |>
    distinct(mirbase_id, .keep_all = TRUE) |> dplyr::select(-species)
}

# Mouse: per-family TargetScanMouse files plus the tables written by
# extract_targets_mesc.R. A family may combine several TargetScan families and
# table labels (the benchmarked mESC set combines three); scores are summed over
# TargetScan families per gene, as in the original mESC benchmark.
load_mouse_evidence <- function(fam, universe, utr) {
  ev <- empty_evidence()
  ts_paths <- file.path(MOUSE_DB$targetscan_dir, fam$targetscan_files)
  found_files <- file.exists(ts_paths)
  if (!all(found_files)) cat("*** TargetScan file(s) not found (skipped):", paste(basename(ts_paths[!found_files]), collapse = ", "), "\n")
  if (any(found_files)) {
    per_family <- map(ts_paths[found_files], function(p) {
      d <- readr::read_tsv(p, show_col_types = FALSE)
      found <- targetscan_column_names(names(d), basename(p))
      gene <- transcript_to_gene(d$`Representative transcript`, utr, universe)
      gene[is.na(gene)] <- symbol_to_gene(d$`Target gene`[is.na(gene)], utr, universe)
      standardise_targetscan(d, found, gene)
    })
    # several TargetScan families (the combined mESC set): scores and site counts summed per
    # gene, as in the original mESC benchmark; PCT = the best family's
    # (na.rm: a family without that information, e.g. no conservation for miR-292a, adds nothing)
    sum_or_na <- function(x) if (all(is.na(x))) NA_real_ else sum(x, na.rm = TRUE)
    ev$TargetScan <- bind_rows(per_family) |> group_by(gene) |>
      summarise(score = sum(score), context = sum_or_na(context), pct = suppressWarnings(max(pct, na.rm = TRUE)),
                n8 = sum_or_na(n8), n7m8 = sum_or_na(n7m8), n7a1 = sum_or_na(n7a1), n_conserved = sum_or_na(n_conserved),
                .groups = "drop") |>
      mutate(pct = ifelse(is.finite(pct), pct, NA_real_))
  }
  td <- MOUSE_DB$tables_dir
  if (!missing_db(file.path(td, "targets_tarbase.tsv"), "TarBase"))
    ev$TarBase <- readr::read_tsv(file.path(td, "targets_tarbase.tsv"), show_col_types = FALSE) |>
      filter(query_mirna %in% fam$table_labels, gene_id %in% universe) |>
      group_by(gene = gene_id) |> summarise(support = suppressWarnings(max(n_records)), .groups = "drop")
  if (!missing_db(file.path(td, "targets_mirtarbase.tsv"), "miRTarBase")) {
    mt <- readr::read_tsv(file.path(td, "targets_mirtarbase.tsv"), show_col_types = FALSE) |>
      filter(query_mirna %in% fam$table_labels) |>
      separate_rows(support_type, sep = ";\\s*") |>
      filter(support_type %in% c("Functional MTI", "Functional MTI (Weak)")) |>
      transmute(gene = symbol_to_gene(gene_symbol, utr, universe), support = support_type, experiments) |> filter(!is.na(gene))
    ev$miRTarBase <- mirtarbase_evidence(mt)
  }
  if (!missing_db(file.path(td, "targets_mirdb.tsv"), "miRDB"))
    ev$miRDB <- readr::read_tsv(file.path(td, "targets_mirdb.tsv"), show_col_types = FALSE) |>
      filter(query_mirna %in% fam$table_labels, !is.na(gene_symbol)) |>
      mutate(gene = symbol_to_gene(gene_symbol, utr, universe)) |> filter(!is.na(gene)) |>
      group_by(gene) |> summarise(score = suppressWarnings(max(score)), .groups = "drop")
  ev
}
