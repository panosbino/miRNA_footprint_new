# ==============================================================================
# HEK293 benchmark, step 0: build every target set in the plan for hsa-miR-124-3p.
#
# All sets are Ensembl gene IDs restricted to one shared gene universe U:
#   genes in scran_normalized_linear.rds (mean > 0) AND present in miReact's
#   human 3'UTR table (tools/miReact/seqs/hs.utr3.seqs.rds).
# Every list-based method (our 5 variants, miReact, miTEA-HiRes) scores the same genes.
#
# Inputs:
#   resources/human/Targets__combined_124-3p_124-3p.2_506-3p.rds      TargetScan
#   resources/human/hsa-miR-124/Targets_miR-124-3p_tarBase.tsv         TarBase v9 (no header)
#   resources/human/hsa-miR-124/Targets_miR-124-3p_miRTarBase.csv      miRTarBase (no header)
#   resources/human/hsa-miR-124/Targets_miR-124-3p_mirDB.txt           NOTE: this is a mirDIP export
#   resources/mouse/miRDB_v6.0_prediction_result.txt.gz                 miRDB v6 (all species; hsa rows used)
# Needs Bioconductor org.Hs.eg.db (miRDB RefSeq -> Ensembl).
#
# Outputs (analysis/HEK_SS3/benchmark/):
#   target_sets.rds, target_sets_summary.tsv, target_sets_symbols.tsv
#
# Environment (all optional):
#   MIRNA_BASE_DIR
#   COEXP_EDGE_PCT        % of candidate pairs kept as edges            (default 0.1)
#   COEXP_MIN_EDGES       minimum edges for a gene to be kept           (default 5)
#   COEXP_NETWORK_CELLS   "induced" (default) or "all"
#   COEXP_MAX_ZERO_FRAC   candidates must be zero in < this fraction of cells (default 0.75)
# ==============================================================================

suppressPackageStartupMessages({ library(tidyverse); library(org.Hs.eg.db) })

BASE_DIR    <- Sys.getenv("MIRNA_BASE_DIR", "/cfs/klemming/projects/supr/naiss2024-6-235/miRNA_footprint_new")
MIREACT_DIR <- file.path(BASE_DIR, "tools/miReact")
PROC_DIR    <- file.path(BASE_DIR, "datasets/HEK_SS3/processed")
HUMAN_DIR   <- file.path(BASE_DIR, "resources/human")
OUT_DIR     <- file.path(BASE_DIR, "analysis/HEK_SS3/benchmark")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
setwd(BASE_DIR)

EDGE_PCT      <- as.numeric(Sys.getenv("COEXP_EDGE_PCT", "0.1"))
MIN_EDGES     <- as.integer(Sys.getenv("COEXP_MIN_EDGES", "5"))
NETWORK_CELLS <- Sys.getenv("COEXP_NETWORK_CELLS", "induced")
MAX_ZERO_FRAC <- as.numeric(Sys.getenv("COEXP_MAX_ZERO_FRAC", "0.75"))
stopifnot(NETWORK_CELLS %in% c("induced", "all"))

MIRNA       <- "hsa-miR-124-3p"
TS_TOP_N    <- c(1000, 500, 200)
MIRDB_TOP   <- c(50, 20, 10, 5)
TARBASE_MIN_RECORDS <- 5

# ---- Data and gene universe ---------------------------------------------------------
counts <- as.matrix(readRDS(file.path(PROC_DIR, "scran_normalized_linear.rds")))
gfp    <- read.csv(file.path(PROC_DIR, "gfp_correlation_data.csv"), stringsAsFactors = FALSE)
stopifnot(setequal(colnames(counts), gfp$cell_id))

utr <- readRDS(file.path(MIREACT_DIR, "seqs/hs.utr3.seqs.rds"))
U <- intersect(rownames(counts)[rowMeans(counts) > 0], utr$gid)
cat(sprintf("scran matrix: %d genes x %d cells; universe U (with 3'UTR) = %d genes\n",
            nrow(counts), ncol(counts), length(U)))

sym_map <- utr |> distinct(gid, gsym) |>
  filter(gid %in% U, !is.na(gsym), gsym != "") |>
  mutate(key = toupper(gsym))
to_gid <- function(symbols) unique(sym_map$gid[sym_map$key %in% toupper(symbols)])
report_mapping <- function(db, symbols) {
  symbols <- unique(na.omit(symbols)); n_ok <- sum(toupper(symbols) %in% sym_map$key)
  cat(sprintf("%s: %d unique symbols, %d map to a gene in U (%d not)\n", db, length(symbols), n_ok, length(symbols) - n_ok))
}

# ---- TargetScan ------------------------------------------------------------------------
# HEK convention (as in every earlier HEK script): rows ordered by cumulative
# weighted context++ score, first occurrence per gene kept (= its best score over
# miR-124-3p.1 / 124-3p.2 / 506-3p). Top N = the N best genes inside U.
ts <- readRDS(file.path(HUMAN_DIR, "Targets__combined_124-3p_124-3p.2_506-3p.rds"))
ts_ranked <- ts |> arrange(Cumulative.weighted.context...score, ensembl_gene_id) |>
  pull(ensembl_gene_id) |> unique()
ts_ranked <- ts_ranked[ts_ranked %in% U]
sets <- list(TargetScan_all = ts_ranked)
for (n in TS_TOP_N) sets[[paste0("TargetScan_top", n)]] <- head(ts_ranked, n)

# ---- TarBase v9 --------------------------------------------------------------------------
# One row per site per experiment; support = number of rows for the gene.
tb <- read_tsv(file.path(HUMAN_DIR, "hsa-miR-124/Targets_miR-124-3p_tarBase.tsv"), col_names = FALSE,
               show_col_types = FALSE)
stopifnot(ncol(tb) == 22, all(tb$X2 == MIRNA))
cat(sprintf("TarBase regulation values: %s\n", paste(unique(tb$X14), collapse = ", ")))
tb_gene <- tb |> mutate(gid = sub("\\..*$", "", X5)) |>
  filter(gid %in% U) |>
  count(gid, name = "support")
sets$TarBase_all <- tb_gene$gid
sets[[sprintf("TarBase_support%d", TARBASE_MIN_RECORDS)]] <- tb_gene$gid[tb_gene$support >= TARBASE_MIN_RECORDS]
sets$TarBase_TSsite <- intersect(tb_gene$gid, ts_ranked)

# ---- miRTarBase ------------------------------------------------------------------------------
# all       "Functional MTI" or "Functional MTI (Weak)" (Non-Functional entries excluded)
# high      "Functional MTI" (strong evidence)
# very high strong evidence including a reporter assay AND Western blot or qRT-PCR
mt <- read_csv(file.path(HUMAN_DIR, "hsa-miR-124/Targets_miR-124-3p_miRTarBase.csv"),
               col_names = c("mti_id", "mirna", "species_mirna", "gene_symbol", "entrez",
                             "species_target", "experiments", "support_type", "pmid"),
               show_col_types = FALSE)
stopifnot(all(mt$mirna == MIRNA))
print(count(mt, support_type))
mt_fun <- filter(mt, support_type %in% c("Functional MTI", "Functional MTI (Weak)"))
strong <- filter(mt, support_type == "Functional MTI")
vh_genes <- strong |> group_by(gene_symbol) |>
  summarise(ok = any(grepl("reporter", experiments, ignore.case = TRUE)) &
                 any(grepl("western|immunoblot|qRT-PCR|qPCR", experiments, ignore.case = TRUE)), .groups = "drop") |>
  filter(ok) |> pull(gene_symbol)
report_mapping("miRTarBase", mt_fun$gene_symbol)
sets$miRTarBase_all      <- to_gid(mt_fun$gene_symbol)
sets$miRTarBase_high     <- to_gid(strong$gene_symbol)
sets$miRTarBase_veryhigh <- to_gid(vh_genes)

# ---- mirDIP (file is named *_mirDB.txt but holds a mirDIP export) ----------------------------
# Score classes: Very High = top 1%, High = top 5%, Medium = top 1/3, Low = rest.
dip <- read_tsv(file.path(HUMAN_DIR, "hsa-miR-124/Targets_miR-124-3p_mirDB.txt"), show_col_types = FALSE)
stopifnot(all(c("Gene Symbol", "MicroRNA", "Score Class") %in% names(dip)), all(dip$MicroRNA == MIRNA))
report_mapping("mirDIP", dip$`Gene Symbol`)
dip_cls <- function(cls) to_gid(dip$`Gene Symbol`[dip$`Score Class` %in% cls])
sets$mirDIP_medium_up <- dip_cls(c("Medium", "High", "Very High"))
sets$mirDIP_high_up   <- dip_cls(c("High", "Very High"))
sets$mirDIP_veryhigh  <- dip_cls("Very High")

# ---- miRDB v6 (hsa rows of the all-species file) --------------------------------------------
md <- read_tsv(file.path(BASE_DIR, "resources/mouse/miRDB_v6.0_prediction_result.txt.gz"),
               col_names = c("mirna", "refseq", "score"), col_types = "ccd") |>
  filter(mirna == MIRNA) |>
  mutate(refseq = sub("\\.\\d+$", "", refseq))
rs_map <- suppressMessages(AnnotationDbi::select(org.Hs.eg.db, keys = unique(md$refseq),
                                                 keytype = "REFSEQ", columns = "ENSEMBL"))
md_gene <- md |> inner_join(rs_map, by = c("refseq" = "REFSEQ"), relationship = "many-to-many") |>
  filter(ENSEMBL %in% U) |>
  group_by(gid = ENSEMBL) |> summarise(score = max(score), .groups = "drop") |>
  arrange(desc(score), gid)
cat(sprintf("miRDB: %d RefSeq predictions -> %d genes in U\n", nrow(md), nrow(md_gene)))
sets$miRDB_all <- md_gene$gid
for (p in MIRDB_TOP) sets[[paste0("miRDB_top", p, "pct")]] <- head(md_gene$gid, ceiling(p / 100 * nrow(md_gene)))

# ---- Co-expression filter -----------------------------------------------------------------------
# Same rule as the mESC benchmark: detection filter (zeros in < MAX_ZERO_FRAC of
# all cells), Spearman correlations among candidates in the network cells, edges
# = top EDGE_PCT % of pairs, keep genes with >= MIN_EDGES edges.
# Network cells: Dox-induced cells by default (as in the HEK co-expression
# scripts). With all cells the induced-vs-uninduced contrast itself, which drives
# GFP, would create the correlations.
zero_frac <- rowMeans(counts == 0)
detected  <- names(zero_frac)[zero_frac < MAX_ZERO_FRAC]
cat(sprintf("\nDetection filter for co-expression candidates: %d of %d genes in U have zeros in < %.0f%% of cells\n",
            length(intersect(detected, U)), length(U), 100 * MAX_ZERO_FRAC))
net_cells <- if (NETWORK_CELLS == "induced") gfp$cell_id[gfp$Dox_Group == "Dox Induced"] else gfp$cell_id

coexp_filter <- function(genes, label) {
  n_in <- length(intersect(genes, U))
  M <- counts[intersect(intersect(genes, U), detected), net_cells, drop = FALSE]
  M <- M[apply(M, 1, var) > 0, , drop = FALSE]
  n <- nrow(M); npairs <- n * (n - 1) / 2
  C <- cor(apply(M, 1, rank, ties.method = "average")); diag(C) <- NA
  ut <- C[upper.tri(C)]
  k <- max(1, ceiling(EDGE_PCT / 100 * npairs))
  thr <- sort(ut, partial = length(ut) - k + 1)[length(ut) - k + 1]
  rm(ut)
  deg <- rowSums(C >= thr, na.rm = TRUE)
  keep <- names(deg)[deg >= MIN_EDGES]
  rm(C); gc(verbose = FALSE)
  cat(sprintf("Co-expression %-17s: %5d of %5d genes pass detection filter; %8.0f pairs, top %.2f%% -> rho >= %.3f, %d edges; %d genes with >= %d edges\n",
              label, n, n_in, npairs, EDGE_PCT, thr, as.integer(sum(deg) / 2), length(keep), MIN_EDGES))
  if (length(keep) < 10) cat(sprintf("  *** WARNING: %s co-expression set has only %d genes ***\n", label, length(keep)))
  list(genes = keep, n_candidates = n, rho_threshold = thr, n_edges = sum(deg) / 2)
}
cat(sprintf("Co-expression networks built on %s cells (n = %d)\n", NETWORK_CELLS, length(net_cells)))
coexp <- list(
  TargetScan        = coexp_filter(sets$TargetScan_all,    "TargetScan"),
  TargetScan_top500 = coexp_filter(sets$TargetScan_top500, "TargetScan_top500"),
  TarBase           = coexp_filter(sets$TarBase_all,       "TarBase"),
  miRTarBase        = coexp_filter(sets$miRTarBase_all,    "miRTarBase"),
  mirDIP            = coexp_filter(sets$mirDIP_medium_up,  "mirDIP")
)
for (db in names(coexp)) sets[[paste0(db, "_coexp")]] <- coexp[[db]]$genes

# ---- Summary / outputs ---------------------------------------------------------------------------
stopifnot(all(map_lgl(sets, \(g) all(g %in% U))))
focal_union <- unique(unlist(sets))
summary_tbl <- tibble(target_set = names(sets), n_genes = lengths(sets)) |> mutate(database = sub("_.*", "", target_set))
cat("\n=== Target sets (genes in U) ===\n"); print(summary_tbl, n = Inf)
if (any(summary_tbl$n_genes < 10)) cat("*** Sets with < 10 genes will be scored but should not be ranked. ***\n")

gid_to_gsym <- utr |> distinct(gid, gsym) |> filter(!is.na(gsym), gsym != "") |> filter(!duplicated(gid))
symbols_tbl <- imap_dfr(sets, \(g, nm) tibble(target_set = nm,
                        gene_symbol = unique(na.omit(gid_to_gsym$gsym[match(g, gid_to_gsym$gid)]))))
write_tsv(summary_tbl, file.path(OUT_DIR, "target_sets_summary.tsv"))
write_tsv(symbols_tbl, file.path(OUT_DIR, "target_sets_symbols.tsv"))
saveRDS(list(sets = sets, universe = U, focal_union = focal_union,
             coexp = map(coexp, \(x) x[c("n_candidates", "rho_threshold", "n_edges")]),
             params = list(edge_pct = EDGE_PCT, min_edges = MIN_EDGES, network_cells = NETWORK_CELLS,
                           max_zero_frac = MAX_ZERO_FRAC, tarbase_min_records = TARBASE_MIN_RECORDS)),
        file.path(OUT_DIR, "target_sets.rds"))
cat(sprintf("\nWrote target_sets.rds, target_sets_summary.tsv, target_sets_symbols.tsv to %s\n", OUT_DIR))
