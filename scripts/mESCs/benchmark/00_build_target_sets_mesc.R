# ==============================================================================
# mESC benchmark, step 0: build every target set in the plan.
#
# All sets are the COMBINED targets of the three mESC families (MIR17, MIR291,
# MIR292), as Ensembl gene IDs, restricted to one shared gene universe U:
#   expressed genes (mean > 0) in final_counts.rds  AND  present in miReact's
#   mouse 3'UTR table (tools/miReact/seqs/mm.utr3.seqs.rds).
# Every list-based method (our 5 variants, miReact, miTEA-HiRes) therefore
# scores exactly the same genes; miReact cannot see genes outside that table.
#
# miRDIP is not built: it is human-only (same situation as miDGD).
#
# Inputs (already in the repo / on Dardel):
#   resources/mouse/Targets_combined_top3_families_mESCs.rds   (TargetScan)
#   resources/mouse/Targets_Negative_Control_mESCs.rds         (reference only)
#   analysis/mESCs/targets/targets_{tarbase,mirtarbase,mirdb}.tsv
#                                         (written by scripts/extract_targets_mesc.R)
#
# Outputs (analysis/mESCs/benchmark/):
#   target_sets.rds          list(sets = named list of Ensembl IDs, universe, focal_union, coexp = network info)
#   target_sets_summary.tsv  one row per set: size and definition
#   target_sets_symbols.tsv  set, gene_symbol   (input for miTEA-HiRes)
#
# Environment variables (all optional):
#   MIRNA_BASE_DIR        project root
#   COEXP_EDGE_PCT        % of candidate pairs kept as edges    (default 0.1)
#   COEXP_MIN_EDGES       minimum edges for a gene to be kept    (default 5)
#   COEXP_NETWORK_CELLS   "control" (default) or "all"
#   COEXP_MAX_ZERO_FRAC   candidates must be zero in < this fraction of cells (default 0.75)
# ==============================================================================

suppressPackageStartupMessages(library(tidyverse))

BASE_DIR    <- Sys.getenv("MIRNA_BASE_DIR", "/cfs/klemming/projects/supr/naiss2024-6-235/miRNA_footprint_new")
MIREACT_DIR <- file.path(BASE_DIR, "tools/miReact")
TARGET_DIR  <- file.path(BASE_DIR, "analysis/mESCs/targets")
OUT_DIR     <- file.path(BASE_DIR, "analysis/mESCs/benchmark")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
setwd(BASE_DIR)

EDGE_PCT      <- as.numeric(Sys.getenv("COEXP_EDGE_PCT", "0.1"))
MIN_EDGES     <- as.integer(Sys.getenv("COEXP_MIN_EDGES", "5"))
NETWORK_CELLS <- Sys.getenv("COEXP_NETWORK_CELLS", "control")
MAX_ZERO_FRAC <- as.numeric(Sys.getenv("COEXP_MAX_ZERO_FRAC", "0.75"))
stopifnot(NETWORK_CELLS %in% c("control", "all"))

FAMILIES    <- c("miR17", "miR291a", "miR292a")  # query_mirna labels in the extracted tables
TS_TOP_N    <- c(1000, 500, 200)
MIRDB_TOP   <- c(50, 20, 10, 5)                   # percent
TARBASE_MIN_RECORDS <- 5

# ---- Data and gene universe ----------------------------------------------------
pheno <- read.csv("datasets/mESCs/phenotype.csv", fileEncoding = "UTF-8-BOM")
if (!"cells" %in% names(pheno)) stop("phenotype.csv: no 'cells' column -- run with LANG=C.UTF-8 (see Pitfalls).")
pheno$Exp[pheno$Exp == "WT"] <- "Control"

counts <- as.matrix(readRDS("datasets/mESCs/final_counts.rds"))
counts <- counts[, colnames(counts) %in% pheno$cells]

utr <- readRDS(file.path(MIREACT_DIR, "seqs/mm.utr3.seqs.rds"))
U <- intersect(rownames(counts)[rowMeans(counts) > 0], utr$gid)
cat(sprintf("Expression: %d genes x %d cells; universe U (expressed, with 3'UTR) = %d genes\n",
            nrow(counts), ncol(counts), length(U)))

# symbol -> Ensembl, case-insensitive, restricted to U (a symbol may map to >1 gene; all are kept)
sym_map <- utr |> distinct(gid, gsym) |>
  filter(gid %in% U, !is.na(gsym), gsym != "") |>
  mutate(key = toupper(gsym))
to_gid <- function(symbols) unique(sym_map$gid[sym_map$key %in% toupper(symbols)])
report_mapping <- function(db, symbols) {
  symbols <- unique(na.omit(symbols))
  n_ok <- sum(toupper(symbols) %in% sym_map$key)
  cat(sprintf("%s: %d unique symbols, %d map to a gene in U (%d not)\n", db, length(symbols), n_ok, length(symbols) - n_ok))
}

# ---- TargetScan ------------------------------------------------------------------
# Same convention as mESC_our_method_separation.R: a gene's score is the SUM of
# its cumulative weighted context++ scores over the three families (more
# negative = stronger). Top N = the N strongest genes INSIDE U, so every set
# really has N genes (the existing top-200 took the first 200 before
# intersecting with the data).
ts <- readRDS("resources/mouse/Targets_combined_top3_families_mESCs.rds") |>
  group_by(ensembl_gene_id) |>
  summarise(score = sum(Cumulative.weighted.context...score), .groups = "drop") |>
  filter(ensembl_gene_id %in% U) |>
  arrange(score, ensembl_gene_id)            # gene ID breaks ties, so the cut is reproducible
ts_ranked <- ts$ensembl_gene_id

sets <- list(TargetScan_all = ts_ranked)
for (n in TS_TOP_N) sets[[paste0("TargetScan_top", n)]] <- head(ts_ranked, n)

# ---- TarBase v9 --------------------------------------------------------------------
# Gene ID column is Ensembl already. Support = number of TarBase records for the
# interaction (independent experiments/samples); a gene's support is the MAX
# over member miRNAs and families (summing would count one CLIP dataset several
# times, once per family member it was assigned to).
tb <- read_tsv(file.path(TARGET_DIR, "targets_tarbase.tsv"), show_col_types = FALSE) |>
  filter(query_mirna %in% FAMILIES, gene_id %in% U)
cat(sprintf("TarBase regulation values kept: %s\n", paste(unique(tb$regulation), collapse = ", ")))
tb_gene <- tb |> group_by(gene_id) |> summarise(support = max(n_records), .groups = "drop")

sets$TarBase_all <- tb_gene$gene_id
sets[[sprintf("TarBase_support%d", TARBASE_MIN_RECORDS)]] <- tb_gene$gene_id[tb_gene$support >= TARBASE_MIN_RECORDS]
sets$TarBase_TSsite <- intersect(tb_gene$gene_id, ts_ranked)     # CLIP-supported AND a predicted seed site

# ---- miRTarBase ----------------------------------------------------------------------
# miRTarBase has no official confidence tiers; these follow its own evidence labels.
#   all        any entry ("Functional MTI" or "Functional MTI (Weak)")
#   high       "Functional MTI" = strong evidence (reporter assay, Western blot, qPCR)
#   very high  strong evidence that includes a reporter assay (direct binding)
#              AND an independent Western blot or qRT-PCR readout
mt <- read_tsv(file.path(TARGET_DIR, "targets_mirtarbase.tsv"), show_col_types = FALSE) |>
  filter(query_mirna %in% FAMILIES) |>
  mutate(strong   = str_detect(support_type, "Functional MTI(;|$)"),
         reporter = str_detect(experiments, regex("reporter", ignore_case = TRUE)),
         readout  = str_detect(experiments, regex("western|immunoblot|qRT-PCR|qPCR", ignore_case = TRUE)))
report_mapping("miRTarBase", mt$gene_symbol)
sets$miRTarBase_all      <- to_gid(mt$gene_symbol)
sets$miRTarBase_high     <- to_gid(mt$gene_symbol[mt$strong])
sets$miRTarBase_veryhigh <- to_gid(mt$gene_symbol[mt$strong & mt$reporter & mt$readout])
print(count(filter(mt, strong), query_mirna, name = "strong_evidence_rows"))

# ---- miRDB ---------------------------------------------------------------------------
# Gene score = max miRDB target score over member miRNAs and families.
# Top X% = the best-scoring X% of the combined list (within U).
md <- read_tsv(file.path(TARGET_DIR, "targets_mirdb.tsv"), show_col_types = FALSE) |>
  filter(query_mirna %in% FAMILIES, !is.na(gene_symbol))
report_mapping("miRDB", md$gene_symbol)
md_gene <- md |>
  mutate(key = toupper(gene_symbol)) |>
  inner_join(dplyr::select(sym_map, gid, key), by = "key", relationship = "many-to-many") |>
  group_by(gid) |> summarise(score = max(score), .groups = "drop") |>
  arrange(desc(score), gid)
sets$miRDB_all <- md_gene$gid
for (p in MIRDB_TOP) sets[[paste0("miRDB_top", p, "pct")]] <- head(md_gene$gid, ceiling(p / 100 * nrow(md_gene)))

# ---- Co-expression filter ----------------------------------------------------------------
# Plan definition: all pairwise Spearman correlations among a database's
# candidates (for TargetScan also the top 500 as a second starting pool); an edge for every pair in the top EDGE_PCT % of correlations;
# keep genes with >= MIN_EDGES edges.
# Network cells: Control only by default. KO cells lack mature miRNAs, so with
# all cells the KO-vs-Control shift itself creates correlation among
# de-repressed targets -- the filter would then select on the very contrast the
# benchmark scores.
# Detection filter (same rule as the HEK scran matrix: zeros in < 75% of cells,
# computed over ALL cells, so it uses no KO/Control label). Without it, the top
# correlations in this 60%-zero CPM matrix are dominated by sparse genes that are
# merely detected in the same cells, and the filter keeps near-unexpressed genes.
zero_frac <- rowMeans(counts == 0)
detected  <- names(zero_frac)[zero_frac < MAX_ZERO_FRAC]
cat(sprintf("\nDetection filter for co-expression candidates: %d of %d genes in U have zeros in < %.0f%% of cells\n",
            length(intersect(detected, U)), length(U), 100 * MAX_ZERO_FRAC))
net_cells <- if (NETWORK_CELLS == "control") pheno$cells[pheno$Exp == "Control"] else pheno$cells
net_cells <- intersect(net_cells, colnames(counts))

coexp_filter <- function(genes, label) {
  n_in <- length(intersect(genes, U))
  M <- counts[intersect(intersect(genes, U), detected), net_cells, drop = FALSE]
  M <- M[apply(M, 1, var) > 0, , drop = FALSE]                 # constant genes have no correlation
  n <- nrow(M); npairs <- n * (n - 1) / 2
  R <- apply(M, 1, rank, ties.method = "average")               # cells x genes; average ranks = Spearman
  C <- cor(R); rm(R)
  diag(C) <- NA
  ut <- C[upper.tri(C)]
  k <- max(1, ceiling(EDGE_PCT / 100 * npairs))
  thr <- sort(ut, partial = length(ut) - k + 1)[length(ut) - k + 1]   # k-th largest correlation
  rm(ut)
  deg <- rowSums(C >= thr, na.rm = TRUE)
  keep <- names(deg)[deg >= MIN_EDGES]
  rm(C); gc(verbose = FALSE)
  cat(sprintf("Co-expression %-12s: %5d of %5d genes pass detection filter\n", label, nrow(M) + 0L, n_in))
  cat(sprintf("Co-expression %-12s: %5d candidates, %9.0f pairs, top %.2f%% -> rho >= %.3f, %d edges; %d genes with >= %d edges\n",
              label, n, npairs, EDGE_PCT, thr, as.integer(sum(deg) / 2), length(keep), MIN_EDGES))
  if (length(keep) < 10) cat(sprintf("  *** WARNING: %s co-expression set has only %d genes ***\n", label, length(keep)))
  list(genes = keep, n_candidates = n, rho_threshold = thr, n_edges = sum(deg) / 2, degree = deg)
}

cat(sprintf("\nCo-expression networks built on %s cells (n = %d)\n", NETWORK_CELLS, length(net_cells)))
coexp <- list(
  TargetScan        = coexp_filter(sets$TargetScan_all,    "TargetScan"),
  TargetScan_top500 = coexp_filter(sets$TargetScan_top500, "TargetScan_top500"),
  TarBase    = coexp_filter(sets$TarBase_all,    "TarBase"),
  miRTarBase = coexp_filter(sets$miRTarBase_all, "miRTarBase")
)
for (db in names(coexp)) sets[[paste0(db, "_coexp")]] <- coexp[[db]]$genes

# ---- Reference: existing TargetScan negative control -------------------------------------
# Targets of miRNA families not expressed in mESCs (from the existing pipeline).
# NOTE: it shares genes with the mESC target lists; that overlap is reported, not removed,
# so results stay comparable with the earlier runs.
negctrl <- unique(readRDS("resources/mouse/Targets_Negative_Control_mESCs.rds")$ensembl_gene_id)
sets$negctrl_TargetScan <- intersect(negctrl, U)

# ---- Summary ----------------------------------------------------------------------------------
stopifnot(all(map_lgl(sets, \(g) all(g %in% U))))
focal_union <- unique(unlist(sets[!grepl("^negctrl", names(sets))]))
cat(sprintf("\nnegctrl_TargetScan: %d genes, of which %d (%.0f%%) are also TargetScan targets of the mESC families\n",
            length(sets$negctrl_TargetScan), sum(sets$negctrl_TargetScan %in% sets$TargetScan_all),
            100 * mean(sets$negctrl_TargetScan %in% sets$TargetScan_all)))

summary_tbl <- tibble(target_set = names(sets), n_genes = lengths(sets)) |>
  mutate(database = sub("_.*", "", target_set))
cat("\n=== Target sets (genes in U) ===\n"); print(summary_tbl, n = Inf)
if (any(summary_tbl$n_genes < 10)) cat("*** Sets with < 10 genes will be scored but should not be ranked. ***\n")

# Symbols for miTEA-HiRes, using the same gid -> symbol rule as export_for_mitea_mesc.R
gid_to_gsym <- utr |> distinct(gid, gsym) |> filter(!is.na(gsym), gsym != "") |> filter(!duplicated(gid))
symbols_tbl <- imap_dfr(sets, \(g, nm) tibble(target_set = nm, gene_symbol = unique(na.omit(gid_to_gsym$gsym[match(g, gid_to_gsym$gid)]))))

write_tsv(summary_tbl, file.path(OUT_DIR, "target_sets_summary.tsv"))
write_tsv(symbols_tbl, file.path(OUT_DIR, "target_sets_symbols.tsv"))
saveRDS(list(sets = sets, universe = U, focal_union = focal_union,
             coexp = map(coexp, \(x) x[c("n_candidates", "rho_threshold", "n_edges")]),
             params = list(edge_pct = EDGE_PCT, min_edges = MIN_EDGES, network_cells = NETWORK_CELLS,
                           max_zero_frac = MAX_ZERO_FRAC,
                           tarbase_min_records = TARBASE_MIN_RECORDS)),
        file.path(OUT_DIR, "target_sets.rds"))
cat(sprintf("\nWrote target_sets.rds, target_sets_summary.tsv, target_sets_symbols.tsv to %s\n", OUT_DIR))
