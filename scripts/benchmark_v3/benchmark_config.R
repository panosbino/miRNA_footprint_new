# ==============================================================================
# Benchmark v3 -- configuration for all three datasets.
#
# Every path and setting lives here. The other scripts take the dataset name
# ("mESC", "HEK" or "PSCSR") as their first argument and look it up below.
# To test with other paths without editing this file, point the environment
# variable BENCHMARK_CONFIG_LOCAL to an R file that overrides some of them.
# ==============================================================================

BASE_DIR    <- Sys.getenv("MIRNA_BASE_DIR", "/cfs/klemming/projects/supr/naiss2024-6-235/miRNA_footprint_new")
MIREACT_DIR <- file.path(BASE_DIR, "tools/miReact")
R_LIBRARY   <- file.path(BASE_DIR, "tools/R_library")        # bayesReact, Regmex
SLURM_ACCOUNT   <- "naiss2026-3-153"
SLURM_PARTITION <- "shared"

# ---- How target lists are built (identical for benchmarked and control miRNAs) ----
TARGETSCAN_TOP_N     <- c(1000, 500, 200)    # strongest predicted targets (cumulative weighted context++)
TARGETSCAN_ALT_TOP_N <- c(500, 200)          # top N by the alternative scores (Aggregate PCT, Total context++)
MIRDB_TOP_PERCENT    <- c(50, 20, 10, 5)     # best-scoring share of the miRDB list
TARBASE_MIN_RECORDS  <- 5                    # independent TarBase records per interaction
COEXP_EDGE_PERCENT   <- 0.1                  # co-expression: top % of gene pairs become edges
COEXP_MIN_EDGES      <- 5                    #   keep genes with at least this many edges
COEXP_MAX_ZERO_FRAC  <- 0.75                 #   candidates must be detected in > 25% of cells
MIN_LIST_SIZE        <- 10                   # smaller lists are scored but not ranked

# ---- Controls ------------------------------------------------------------------------------
N_RANDOM    <- 10                            # random gene lists (and random motifs / miRNAs) per list
# Which genes count as targets of the benchmarked family, and are therefore removed
# from its negative-control lists and from its random-gene pool:
#   "all_databases"  genes listed by ANY database (default; the cleanest null)
#   "same_database"  only genes listed by the database the list comes from (larger
#                    pools, but they still contain targets known to other databases)
CONTROL_EXCLUSION <- "all_databases"
# mESC single-family benchmarks: genes shared with the other benchmarked families
# are removed from each family's lists, so a family is tested on its own targets.
#   "same_database"  a gene is shared if the same database lists it for another family
#                    (default; with "all_databases" TarBase alone would remove ~half the genes)
#   "all_databases"  a gene is shared if any database lists it for another family
SHARED_TARGET_EXCLUSION <- "same_database"
RANDOM_SEED <- 20261023

# ---- Human target databases (HEK and PSCSR) ---------------------------------------------------
HUMAN_DB <- list(
  targetscan_summary  = file.path(BASE_DIR, "resources/human/TargetScan8/Summary_Counts.all_predictions.txt"),
  targetscan_families = file.path(BASE_DIR, "resources/human/TargetScan8/miR_Family_Info.txt"),
  tarbase             = file.path(BASE_DIR, "resources/human/Homo_sapiens_TarBase-v9.tsv"),
  mirtarbase          = file.path(BASE_DIR, "resources/human/hsa_MTI.csv"),
  mirdip              = file.path(BASE_DIR, "resources/human/mirDIP_benchmark_v2.txt"),   # same export as v2; must include the control miRNA
  mirdb               = file.path(BASE_DIR, "resources/mouse/miRDB_v6.0_prediction_result.txt.gz")   # all species
)

# ---- Mouse target databases (mESC) --------------------------------------------------------------
MOUSE_DB <- list(
  targetscan_dir = file.path(BASE_DIR, "resources/mouse/Targetscan_files_mouse/mouse"),   # per-family TargetScanMouse 8.0 files
  tables_dir     = file.path(BASE_DIR, "analysis/mESCs/targets"),     # targets_{tarbase,mirtarbase,mirdb}.tsv (extract_targets_mesc.R)
  mature_fasta   = file.path(BASE_DIR, "resources/mouse/mature.fa")
)

# ---- Datasets ----------------------------------------------------------------------------------------
DATASETS <- list(
  mESC = list(
    species          = "mm",
    out_dir          = "analysis/mESCs/benchmark_v3",
    expression       = "datasets/mESCs/final_counts.rds",            # CPM (no raw counts available)
    mitea_counts_dir = "analysis/mESCs/mitea_input/counts_data",
    # The three seed families highly expressed in mESCs (all lost in the Drosha KO).
    # "combined": one list per database for the three together, activities summed for
    # motif methods. The single families are benchmarked on their own targets: genes
    # shared with the other two families are removed (SHARED_TARGET_EXCLUSION).
    benchmark_families = list(
      `miR-17/291/292 (combined)` = list(
        targetscan_files = c("TargetScan8.0__miR-17-5p_20-5p_93-5p_106-5p.predicted_targets.txt",
                             "TargetScan8.0__miR-291-3p_294-3p_295-3p_302-3p.predicted_targets.txt",
                             "TargetScan8.0__miR-292a-3p_467a-5p.predicted_targets.txt"),
        table_labels     = c("miR17", "miR291a", "miR292a"),        # query_mirna in the extracted tables
        seed_mirnas      = c("mmu-miR-17-5p", "mmu-miR-291a-3p", "mmu-miR-292a-3p")),
      `miR-17` = list(
        targetscan_files = "TargetScan8.0__miR-17-5p_20-5p_93-5p_106-5p.predicted_targets.txt",
        table_labels = "miR17", seed_mirnas = "mmu-miR-17-5p"),
      `miR-291a` = list(
        targetscan_files = "TargetScan8.0__miR-291-3p_294-3p_295-3p_302-3p.predicted_targets.txt",
        table_labels = "miR291a", seed_mirnas = "mmu-miR-291a-3p"),
      `miR-292a` = list(
        targetscan_files = "TargetScan8.0__miR-292a-3p_467a-5p.predicted_targets.txt",
        table_labels = "miR292a", seed_mirnas = "mmu-miR-292a-3p")),
    # Negative control: miR-122-5p, liver-specific and not expressed in mESCs. Its
    # TarBase/miRTarBase/miRDB rows come from extract_targets_mesc.R (family "miR122").
    control_family = list(
      name             = "miR-122-5p",
      targetscan_files = "TargetScan8.0__miR-122-5p.predicted_targets.txt",
      table_labels     = "miR122",
      seed_mirnas      = "mmu-miR-122-5p")
  ),
  HEK = list(
    species          = "hs",
    out_dir          = "analysis/HEK_SS3/benchmark_v3",
    expression       = "datasets/HEK_SS3/processed/scran_normalized_linear.rds",
    mitea_counts_dir = "analysis/HEK_SS3/comparisons/mitea_input/counts_data",
    benchmark_mirnas = "hsa-miR-124-3p",                 # its seed family is benchmarked
    control_mirna    = "hsa-miR-122-5p"
  ),
  PSCSR = list(
    species          = "hs",
    out_dir          = "analysis/PSCSR/benchmark_v3",
    expression       = "analysis/PSCSR/benchmark_v3/data/scran_normalized_linear.rds",   # written by pscsr_0
    mitea_counts_dir = "analysis/PSCSR/benchmark_v3/data/mitea_input/counts_data",
    control_mirna    = "hsa-miR-122-5p"
    # benchmarked families: the top PSCSR_TOP_N per cell line, chosen by pscsr_0
  )
)

# ---- PSCSR data preparation (pscsr_0_prepare_data.R) -------------------------------------------------
PSCSR_INPUT <- list(
  dir          = file.path(BASE_DIR, "datasets/SC_PSCSR_v2"),
  mrna_counts  = "hs_PCSCRv2_mrna_counts.tsv",            # Ensembl x cells, UMIs
  mirna_counts = "hs_PCSCRv2_mirna_counts_mirbase.tsv",   # mature miRBase IDs x cells, UMIs
  mirna_anno   = "hs_PCSCRv2_mirna_anno_mirbase.tsv",     # only used to cross-check seed sites
  cell_anno    = "hs_PCSCRv2_cell_anno.tsv",
  cell_col     = "barcode",
  line_col     = "CL_short")
PSCSR_CELL_LINES    <- c("A549", "K562", "HeLa", "293T")
PSCSR_TOP_N         <- 5                  # most expressed miRNA families per cell line
PSCSR_MIRNA_NORM    <- "median_ratio"     # or "cpm"
PSCSR_PROTEIN_CODING <- file.path(BASE_DIR, "datasets/HEK_SS3/processed/biomart_protein_coding_ids.rds")

if (nzchar(Sys.getenv("BENCHMARK_CONFIG_LOCAL"))) source(Sys.getenv("BENCHMARK_CONFIG_LOCAL"))

# ---- Helpers to resolve a dataset -------------------------------------------------------------------------
get_dataset <- function(name) {
  if (!name %in% names(DATASETS)) stop("Unknown dataset '", name, "'. Use one of: ", paste(names(DATASETS), collapse = ", "))
  d <- DATASETS[[name]]
  d$name <- name
  for (k in c("out_dir", "expression", "mitea_counts_dir")) d[[k]] <- file.path(BASE_DIR, d[[k]])
  dir.create(d$out_dir, recursive = TRUE, showWarnings = FALSE)
  d
}
dataset_from_args <- function() {
  a <- commandArgs(trailingOnly = TRUE)
  if (!length(a)) stop("Give the dataset name as the first argument: mESC, HEK or PSCSR")
  get_dataset(a[1])
}
