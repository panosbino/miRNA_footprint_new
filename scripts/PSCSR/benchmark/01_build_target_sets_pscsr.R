# ==============================================================================
# PSCSR benchmark, step 1: target sets for every unit (seed family), same set
# types and rules as the HEK benchmark (00_build_target_sets_hek.R).
#
# Database lookups use ALL human miRNAs of the unit (units.tsv "members"),
# because family members share targets; per gene, the best evidence over
# members is kept (TargetScan: min cumulative weighted context++; TarBase:
# max records per member; miRDB: max score; mirDIP: best class).
# TargetScan rows are matched to a unit by SEED: in TargetScan 8.0's Summary_Counts
# the family column holds the seed (nt 2-8, RNA letters), not the family name.
#
# Gene universe U: genes in the prepared scran matrix (mean > 0) with a miReact
# human 3'UTR entry.
#
# Co-expression network cells (COEXP_NETWORK_CELLS):
#   "within" (default) -- Spearman correlation pooled WITHIN cell types (ranks
#            z-scored inside each type, then correlated). Across all cells the
#            four cell lines dominate every correlation, so the filter would
#            pick cell-line marker modules, and those also drive the across-cell
#            evaluation.
#   "all"    -- plain Spearman over all cells.
#
# Outputs (analysis/PSCSR/benchmark/): target_sets.rds, target_sets_summary.tsv,
#   target_sets_symbols.tsv. Set keys are "<unit>|<set type>".
# ==============================================================================

suppressPackageStartupMessages({ library(tidyverse); library(org.Hs.eg.db) })
source(file.path(Sys.getenv("MIRNA_BASE_DIR", "/cfs/klemming/projects/supr/naiss2024-6-235/miRNA_footprint_new"),
                 "scripts/PSCSR/benchmark/config_pscsr.R"))

source(file.path(Sys.getenv("MIRNA_BASE_DIR", "~/Desktop/Projects/miRNA_footprint_new"),
                 "scripts/PSCSR/benchmark/config_pscsr.R"))

EDGE_PCT      <- as.numeric(Sys.getenv("COEXP_EDGE_PCT", "0.1"))
MIN_EDGES     <- as.integer(Sys.getenv("COEXP_MIN_EDGES", "5"))
NETWORK_CELLS <- Sys.getenv("COEXP_NETWORK_CELLS", "within")
stopifnot(NETWORK_CELLS %in% c("within", "all"))
TS_TOP_N  <- c(1000, 500, 200)
MIRDB_TOP <- c(50, 20, 10, 5)
TARBASE_MIN_RECORDS <- 5

counts <- readRDS(file.path(DATA_DIR, "scran_normalized_linear.rds"))
cells  <- read_tsv(file.path(DATA_DIR, "cells.tsv"), show_col_types = FALSE)
units  <- read_tsv(file.path(DATA_DIR, "units.tsv"), show_col_types = FALSE)
stopifnot(identical(colnames(counts), cells$cell))
members <- setNames(strsplit(units$members, ";"), units$unit)

utr <- readRDS(file.path(MIREACT_DIR, "seqs/hs.utr3.seqs.rds"))
U <- intersect(rownames(counts)[rowMeans(counts) > 0], utr$gid)
cat(sprintf("Universe U: %d genes; %d units\n", length(U), nrow(units)))
sym_map <- utr |> distinct(gid, gsym) |> filter(gid %in% U, !is.na(gsym), gsym != "") |> mutate(key = toupper(gsym))
tid_map <- utr |> distinct(tid, gid) |> mutate(tid = sub("\\.\\d+$", "", tid)) |> filter(gid %in% U)

# A database whose file is missing is SKIPPED (its sets stay empty) so the rest can
# run, e.g. before the mirDIP export exists; this is printed loudly.
have <- function(path, db) { ok <- file.exists(path); if (!ok) cat(sprintf("*** %s SKIPPED: file not found: %s\n", db, path)); ok }
read_header <- function(path, ...) names(data.table::fread(path, nrows = 0, ...))
fread_cols  <- function(path, cols, ...) data.table::fread(path, select = cols, data.table = FALSE, ...)

# ---- TargetScan (bulk Summary_Counts, human) -------------------------------------------------
ts_ranked <- setNames(rep(list(character(0)), nrow(units)), units$unit)
if (have(TS_SUMMARY_FILE, "TargetScan")) {
h <- read_header(TS_SUMMARY_FILE)
c_tx <- pick_col(setNames(as.list(h), h), c("Transcript ID"), "TS transcript")
c_sy <- pick_col(setNames(as.list(h), h), c("Gene Symbol"), "TS gene symbol")
c_fm <- pick_col(setNames(as.list(h), h), c("miRNA family", "miR Family", "miR family"), "TS family")
c_sp <- pick_col(setNames(as.list(h), h), c("Species ID"), "TS species")
c_sc <- pick_col(setNames(as.list(h), h), c("Cumulative weighted context++ score"), "TS score")
ts_raw <- fread_cols(TS_SUMMARY_FILE, c(c_tx, c_sy, c_fm, c_sp, c_sc))
names(ts_raw) <- c("tx", "sym", "family", "species", "score")
ts_raw <- ts_raw |> filter(species == 9606) |> mutate(score = suppressWarnings(as.numeric(score)), tx = sub("\\.\\d+$", "", tx))
ts_raw$gid <- tid_map$gid[match(ts_raw$tx, tid_map$tid)]
no_tx <- is.na(ts_raw$gid)
ts_raw$gid[no_tx] <- sym_map$gid[match(toupper(ts_raw$sym[no_tx]), sym_map$key)]
ts_raw$seed <- toupper(gsub("U", "T", ts_raw$family))           # family column = seed in RNA letters
miss_seed <- setdiff(units$seed, ts_raw$seed)
if (length(miss_seed)) cat("*** Unit seeds with no TargetScan rows:", paste(units$unit[units$seed %in% miss_seed], collapse = ", "), "\n")
cat(sprintf("TargetScan: %d human rows, %d genes in U\n", nrow(ts_raw), n_distinct(na.omit(ts_raw$gid))))
ts_ranked <- map(setNames(units$seed, units$unit), function(sd) {
  ts_raw |> filter(seed == sd, !is.na(gid), !is.na(score)) |>
    group_by(gid) |> summarise(score = suppressWarnings(min(score)), .groups = "drop") |>
    arrange(score, gid) |> pull(gid)
})
}

# ---- TarBase v9 (human) ---------------------------------------------------------------------------
tb <- tibble(mirna = character(), gid = character(), records = integer())
if (have(TARBASE_FILE, "TarBase")) {
h <- read_header(TARBASE_FILE)
hl <- setNames(as.list(h), h)
tb <- fread_cols(TARBASE_FILE, c(pick_col(hl, c("mirna_name", "mirna"), "TarBase miRNA"),
                                 pick_col(hl, c("gene_id"), "TarBase gene_id")))
names(tb) <- c("mirna", "gene_id")
tb <- tb |> filter(mirna %in% unlist(members)) |> mutate(gid = sub("\\..*$", "", gene_id)) |> filter(gid %in% U) |>
  count(mirna, gid, name = "records")
}

# ---- miRTarBase (human) -----------------------------------------------------------------------------
mt <- tibble(mirna = character(), gene = character(), species = character(), exp = character(), support = character())
if (have(MIRTARBASE_FILE, "miRTarBase")) {
mt <- read_csv(MIRTARBASE_FILE, show_col_types = FALSE)
mt <- tibble(mirna   = mt[[pick_col(mt, c("miRNA"), "miRTarBase miRNA")]],
             gene    = mt[[pick_col(mt, c("Target Gene"), "miRTarBase target")]],
             species = mt[[pick_col(mt, c("Species (Target Gene)"), "miRTarBase species")]],
             exp     = mt[[pick_col(mt, c("Experiments"), "miRTarBase experiments")]],
             support = mt[[pick_col(mt, c("Support Type"), "miRTarBase support")]]) |>
  filter(species == "hsa", mirna %in% unlist(members))
}

# ---- mirDIP (export for the selected miRNAs) --------------------------------------------------------
dip <- tibble(`Gene Symbol` = character(), MicroRNA = character(), `Score Class` = character()); missing_dip <- character(0)
if (have(MIRDIP_FILE, "mirDIP")) {
dip <- read_tsv(MIRDIP_FILE, show_col_types = FALSE)
stopifnot(all(c("Gene Symbol", "MicroRNA", "Score Class") %in% names(dip)))
missing_dip <- setdiff(unlist(members), dip$MicroRNA)
}

# ---- miRDB v6 (hsa rows) ------------------------------------------------------------------------------
md <- read_tsv(MIRDB_FILE, col_names = c("mirna", "refseq", "score"), col_types = "ccd") |>
  filter(mirna %in% unlist(members)) |> mutate(refseq = sub("\\.\\d+$", "", refseq))
rs_map <- suppressMessages(AnnotationDbi::select(org.Hs.eg.db, keys = unique(md$refseq), keytype = "REFSEQ", columns = "ENSEMBL"))
md <- md |> inner_join(rs_map, by = c("refseq" = "REFSEQ"), relationship = "many-to-many") |> filter(ENSEMBL %in% U)

to_gid <- function(symbols) unique(sym_map$gid[sym_map$key %in% toupper(symbols)])

# ---- Co-expression filter (rule as HEK/mESC) --------------------------------------------------------------
detected <- rownames(counts)[rowMeans(counts == 0) < 0.75]
zrank_within <- function(M) {                       # cells x genes, ranks z-scored inside each cell type
  Z <- matrix(0, ncol(M), nrow(M), dimnames = list(colnames(M), rownames(M)))
  for (ct in unique(cells$cell_type)) {
    i <- which(cells$cell_type == ct)
    R <- apply(M[, i, drop = FALSE], 1, rank)
    Z[i, ] <- scale(R)
  }
  Z[!is.finite(Z)] <- 0                              # genes constant within a type contribute nothing there
  Z
}
coexp_filter <- function(genes, label) {
  M <- counts[intersect(intersect(genes, U), detected), , drop = FALSE]
  M <- M[apply(M, 1, var) > 0, , drop = FALSE]
  n <- nrow(M)
  if (n < 10) { cat(sprintf("  %-34s %4d candidates -- skipped\n", label, n)); return(character(0)) }
  C <- if (NETWORK_CELLS == "within") cor(zrank_within(M)) else cor(apply(M, 1, rank))
  diag(C) <- NA
  ut <- C[upper.tri(C)]
  k <- max(1, ceiling(EDGE_PCT / 100 * length(ut)))
  thr <- sort(ut, partial = length(ut) - k + 1)[length(ut) - k + 1]
  deg <- rowSums(C >= thr, na.rm = TRUE)
  keep <- names(deg)[deg >= MIN_EDGES]
  cat(sprintf("  %-34s %4d candidates, rho >= %.3f, %d edges -> %d genes\n", label, n, thr, k, length(keep)))
  keep
}

# ---- Build per unit ---------------------------------------------------------------------------------------
cat(sprintf("\nCo-expression networks: %s cells\n", NETWORK_CELLS))
sets <- list(); background <- list()
for (u in units$unit) {
  mem <- members[[u]]
  S <- list()
  S$TargetScan_all <- ts_ranked[[u]]
  for (n in TS_TOP_N) S[[paste0("TargetScan_top", n)]] <- head(ts_ranked[[u]], n)

  tbu <- tb |> filter(mirna %in% mem) |> group_by(gid) |> summarise(support = suppressWarnings(max(records)), .groups = "drop")
  S$TarBase_all <- tbu$gid
  S[[sprintf("TarBase_support%d", TARBASE_MIN_RECORDS)]] <- tbu$gid[tbu$support >= TARBASE_MIN_RECORDS]
  S$TarBase_TSsite <- intersect(tbu$gid, ts_ranked[[u]])

  mtu <- mt |> filter(mirna %in% mem)
  strong <- mtu |> filter(support == "Functional MTI")
  vh <- strong |> group_by(gene) |>
    summarise(ok = any(grepl("reporter", exp, ignore.case = TRUE)) &
                   any(grepl("western|immunoblot|qRT-PCR|qPCR", exp, ignore.case = TRUE)), .groups = "drop") |>
    filter(ok) |> pull(gene)
  S$miRTarBase_all      <- to_gid(mtu$gene[mtu$support %in% c("Functional MTI", "Functional MTI (Weak)")])
  S$miRTarBase_high     <- to_gid(strong$gene)
  S$miRTarBase_veryhigh <- to_gid(vh)

  dipu <- dip |> filter(MicroRNA %in% mem)
  S$mirDIP_medium_up <- to_gid(dipu$`Gene Symbol`[dipu$`Score Class` %in% c("Medium", "High", "Very High")])
  S$mirDIP_high_up   <- to_gid(dipu$`Gene Symbol`[dipu$`Score Class` %in% c("High", "Very High")])
  S$mirDIP_veryhigh  <- to_gid(dipu$`Gene Symbol`[dipu$`Score Class` == "Very High"])

  mdu <- md |> filter(mirna %in% mem) |> group_by(gid = ENSEMBL) |> summarise(score = max(score), .groups = "drop") |>
    arrange(desc(score), gid)
  S$miRDB_all <- mdu$gid
  for (p in MIRDB_TOP) S[[paste0("miRDB_top", p, "pct")]] <- head(mdu$gid, ceiling(p / 100 * nrow(mdu)))

  cat(sprintf("\n[%s] co-expression:\n", u))
  S$TargetScan_coexp        <- coexp_filter(S$TargetScan_all,    "TargetScan_all")
  S$TargetScan_top500_coexp <- coexp_filter(S$TargetScan_top500, "TargetScan_top500")
  S$TarBase_coexp           <- coexp_filter(S$TarBase_all,       "TarBase_all")
  S$miRTarBase_coexp        <- coexp_filter(S$miRTarBase_all,    "miRTarBase_all")
  S$mirDIP_coexp            <- coexp_filter(S$mirDIP_medium_up,  "mirDIP_medium_up")

  S <- map(S, \(g) intersect(g, U))
  background[[u]] <- setdiff(U, unlist(S))                     # C2: genes in no target set of this unit
  names(S) <- set_key(u, names(S))
  sets <- c(sets, S)
}
if (length(missing_dip)) cat(sprintf("\n*** Not in the mirDIP file (no mirDIP sets for them): %s\n", paste(missing_dip, collapse = ", ")))

summary_tbl <- tibble(target_set = names(sets), n_genes = lengths(sets)) |> bind_cols(split_key(names(sets)))
cat("\n=== Set sizes (rows = set type, columns = unit) ===\n")
print(summary_tbl |> dplyr::select(set_type, unit, n_genes) |> pivot_wider(names_from = unit, values_from = n_genes) |>
        as.data.frame(), row.names = FALSE)
cat("Sets with < 10 genes are scored but should not be ranked; empty sets are skipped.\n")

gid_to_gsym <- utr |> distinct(gid, gsym) |> filter(!is.na(gsym), gsym != "") |> filter(!duplicated(gid))
symbols_tbl <- imap_dfr(sets, \(g, nm) tibble(target_set = nm,
                        gene_symbol = unique(na.omit(gid_to_gsym$gsym[match(g, gid_to_gsym$gid)]))))
write_tsv(summary_tbl, file.path(OUT_DIR, "target_sets_summary.tsv"))
write_tsv(symbols_tbl, file.path(OUT_DIR, "target_sets_symbols.tsv"))
saveRDS(list(sets = sets, universe = U, background = background,
             params = list(edge_pct = EDGE_PCT, min_edges = MIN_EDGES, network_cells = NETWORK_CELLS,
                           tarbase_min_records = TARBASE_MIN_RECORDS)),
        file.path(OUT_DIR, "target_sets.rds"))
cat(sprintf("Wrote target_sets.rds (%d sets), target_sets_summary.tsv, target_sets_symbols.tsv\n", length(sets)))

