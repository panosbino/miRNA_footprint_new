# ==============================================================================
# Dataset 3 (PSCSR-seq V2, Li et al. 2025 Sci Rep, GSE226714): configuration.
# Cultured A549, K562, HeLa, 293T cells with miRNA and mRNA from the same cell.
#
# EVERYTHING THAT DEPENDS ON THE COUNT TABLES IS IN THIS FILE. When the tables
# arrive, only the "INPUT FILES" block (and, if needed, the two loader
# functions) should have to change. Every other script sources this file.
# ==============================================================================

BASE_DIR    <- Sys.getenv("MIRNA_BASE_DIR", "~/Desktop/Projects/miRNA_footprint_new")
MIREACT_DIR <- file.path(BASE_DIR, "tools/miReact")
OUT_DIR     <- file.path(BASE_DIR, "analysis/PSCSR/benchmark")
DATA_DIR    <- file.path(OUT_DIR, "data")          # prepared matrices (step 0)

# ---- INPUT FILES (datasets/SC_PSCSR_v2, processed data from GSE226714) ----------------
PSCSR_DIR          <- file.path(BASE_DIR, "datasets/SC_PSCSR_v2")
MRNA_COUNTS_FILE   <- file.path(PSCSR_DIR, "hs_PCSCRv2_mrna_counts.tsv")          # Ensembl x cells, UMIs (2.3% non-integer: multi-mapping weights)
MIRNA_COUNTS_FILE  <- file.path(PSCSR_DIR, "hs_PCSCRv2_mirna_counts_mirbase.tsv") # mature miRBase IDs x cells, UMIs
MIRNA_ANNO_FILE    <- file.path(PSCSR_DIR, "hs_PCSCRv2_mirna_anno_mirbase.tsv")   # used only to cross-check motifs (see step 0)
CELL_META_FILE     <- file.path(PSCSR_DIR, "hs_PCSCRv2_cell_anno.tsv")
COUNTS_ORIENTATION <- "features_x_cells"
META_CELL_COL      <- "barcode"
META_TYPE_COL      <- "CL_short"           # A549 / K562 / HeLa / 293T
CELL_TYPES         <- c("A549", "K562", "HeLa", "293T")

# ---- Target databases (human; TODO: paths to the downloaded files) -----------------
# TargetScan 8.0 human bulk files (targetscan.org > Download data):
TS_SUMMARY_FILE   <- file.path(BASE_DIR, "resources/human/TargetScan8/Summary_Counts.all_predictions.txt")
TS_FAMILY_FILE    <- file.path(BASE_DIR, "resources/human/TargetScan8/miR_Family_Info.txt")
TARBASE_FILE      <- file.path(BASE_DIR, "resources/human/Homo_sapiens_TarBase-v9.tsv")   # DIANA TarBase v9, human
MIRTARBASE_FILE   <- file.path(BASE_DIR, "resources/human/hsa_MTI.csv")                  # miRTarBase, human, with header
MIRDIP_FILE       <- file.path(BASE_DIR, "resources/human/mirDIP_PSCSR_selected.txt")    # mirDIP export for the selected miRNAs (step 0 lists them)
MIRDB_FILE        <- file.path(BASE_DIR, "resources/mouse/miRDB_v6.0_prediction_result.txt.gz")  # all species
PROTEIN_CODING_FILE <- file.path(BASE_DIR, "datasets/HEK_SS3/processed/biomart_protein_coding_ids.rds")

# ---- Benchmark parameters -------------------------------------------------------------
TOP_N_PER_TYPE <- 10       # top expressed units per cell line (each line evaluated on its own list)
MAX_ZERO_PROP  <- 0.75     # mRNA gene filter, as for the HEK scran matrix
SCRAN_SEED     <- 123      # as build_annotated_scran_dataset.R
UNIT_LEVEL     <- "family" # "family": miRNAs sharing nt 2-8 are one unit (targets are shared);
                           # measured expression = sum of all detected members
MIRNA_NORM     <- "median_ratio"  # ground-truth normalisation: "median_ratio" (median of ratios over miRNAs
                                  # detected in >= 90% of cells, per cell line) or "cpm". CPM is compositional:
                                  # the miR-17 family is ~1/3 of K562 reads, so it depresses every other miRNA.
SPEC_MAX_TRUTH_RHO <- 0.3  # specificity: other units whose MEASURED expression co-varies with the unit's
                           # own (|within-line rho| >= this, e.g. miR-17~92 cluster members) are not used

# ---- Generic readers ------------------------------------------------------------------------
# A matrix file may be .rds (matrix, dgCMatrix, data.frame) or a delimited text
# file (.csv/.tsv/.txt, optionally .gz) whose first column holds the row names.
read_count_matrix <- function(path) {
  if (!file.exists(path)) stop("Count table not found: ", path, "  (set it in config_pscsr.R)")
  m <- if (grepl("\\.rds$", path, ignore.case = TRUE)) readRDS(path) else {
    sep <- if (grepl("\\.csv(\\.gz)?$", path, ignore.case = TRUE)) "," else "\t"
    # R write.table() files have one header field fewer than data fields; fread then
    # names the row-name column V1 and warns -- expected, so that warning is muffled.
    df <- withCallingHandlers(data.table::fread(path, sep = sep, data.table = FALSE),
                              warning = function(w) if (grepl("column names but the data has", conditionMessage(w))) invokeRestart("muffleWarning"))
    rownames(df) <- df[[1]]; df[[1]] <- NULL; df
  }
  m <- as.matrix(m)
  if (COUNTS_ORIENTATION == "cells_x_features") m <- t(m)
  storage.mode(m) <- "double"
  if (anyDuplicated(rownames(m))) stop(path, ": duplicated feature names")
  if (anyDuplicated(colnames(m))) stop(path, ": duplicated cell names")
  m
}

load_mrna_counts  <- function() read_count_matrix(MRNA_COUNTS_FILE)    # genes x cells, RAW UMI counts
load_mirna_counts <- function() read_count_matrix(MIRNA_COUNTS_FILE)   # miRNAs x cells, RAW UMI counts

load_cell_meta <- function() {
  if (!file.exists(CELL_META_FILE)) stop("Cell metadata not found: ", CELL_META_FILE)
  md <- data.table::fread(CELL_META_FILE, data.table = FALSE)
  if (!all(c(META_CELL_COL, META_TYPE_COL) %in% names(md)))
    stop("Metadata needs columns '", META_CELL_COL, "' and '", META_TYPE_COL, "'; found: ", paste(names(md), collapse = ", "))
  out <- data.frame(cell = as.character(md[[META_CELL_COL]]), cell_type = as.character(md[[META_TYPE_COL]]))
  bad <- setdiff(unique(out$cell_type), CELL_TYPES)
  if (length(bad)) stop("Unexpected cell types in metadata: ", paste(bad, collapse = ", "))
  out
}

# miRNA names -> miRBase mature IDs ("hsa-miR-21-5p"). Adjust here if the table
# uses another convention (e.g. precursor names, which have no arm and cannot be
# mapped to a seed without more information).
standardise_mirna_names <- function(x) {
  x <- trimws(x)
  ifelse(grepl("^hsa-", x), x, paste0("hsa-", x))
}

# ---- Shared helpers ---------------------------------------------------------------------------
pick_col <- function(df, candidates, what = "") {
  hit <- candidates[candidates %in% names(df)]
  if (!length(hit)) stop(sprintf("No column for %s among: %s. Columns are: %s", what,
                                 paste(candidates, collapse = ", "), paste(names(df), collapse = ", ")))
  hit[1]
}
revcomp_dna <- function(s) {
  s <- toupper(gsub("U", "T", s))
  paste(rev(strsplit(chartr("ACGT", "TGCA", s), "")[[1]]), collapse = "")
}
# key of a target set in all flat tables: "<unit>|<set type>"
set_key   <- function(unit, set_type) paste(unit, set_type, sep = "|")
split_key <- function(key) data.frame(unit = sub("\\|.*$", "", key), set_type = sub("^[^|]*\\|", "", key))

# Optional local overrides (e.g. test paths) without editing this file:
if (nzchar(Sys.getenv("PSCSR_CONFIG_LOCAL"))) source(Sys.getenv("PSCSR_CONFIG_LOCAL"))
