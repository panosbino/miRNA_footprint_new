# ==============================================================================
# Shared setup for the TarBase diagnostics (scripts 01-04).
#
# Question: why does our aggregative method do so much worse with the TarBase
# list than ranked methods (miReact, miTEA-HiRes) given the same list?
#
# Sourced by every script in this folder. Loads data once, defines the target
# lists, the scoring variants, covariates, and evaluation shortcuts.
#
# Configuration (environment variables, all optional):
#   MIRNA_BASE_DIR     project root (default: current working directory)
#   MIRNA_TARBASE_RDS  path to miReact's tarbase.rds
#                      (default: <root>/tools/miReact/data/tarbase.rds)
#   DIAG_N_BOOT        bootstrap resamples for confidence intervals (default 200)
#   DIAG_N_RANDOM      random draws for null / subset analyses (default 20)
# ==============================================================================

suppressPackageStartupMessages({ library(tidyverse); library(Matrix) })   # Matrix: raw UMI file is a sparse dgCMatrix

BASE_DIR      <- Sys.getenv("MIRNA_BASE_DIR", getwd())
PROCESSED_DIR <- file.path(BASE_DIR, "datasets/HEK_SS3/processed")
OUT_DIR       <- file.path(BASE_DIR, "analysis/HEK_SS3/tarbase_diagnostics")
TARBASE_RDS   <- Sys.getenv("MIRNA_TARBASE_RDS", file.path(BASE_DIR, "tools/miReact/data/tarbase.rds"))
N_BOOT        <- as.integer(Sys.getenv("DIAG_N_BOOT", "200"))
N_RANDOM      <- as.integer(Sys.getenv("DIAG_N_RANDOM", "20"))
SEED          <- 1312
TARGET_MIRNA  <- "hsa-miR-124-3p"
N_TOP_TS      <- 200
PSEUDOCOUNT   <- 1   # for log ratios; scran-normalized values are on the UMI scale
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

source(file.path(BASE_DIR, "scripts/Utils.R"))   # calculate_activity(), load_gfp_truth(), evaluate_score_vs_gfp()

# ---- Expression ---------------------------------------------------------------
X <- as.matrix(readRDS(file.path(PROCESSED_DIR, "scran_normalized_linear.rds")))   # genes x cells, linear scale
cat(sprintf("scran matrix: %d genes x %d cells\n", nrow(X), ncol(X)))

raw <- readRDS(file.path(PROCESSED_DIR, "combined_counts_filtered_UMIs.rds"))
rownames(raw) <- sub("\\.\\d+$", "", rownames(raw))          # drop Ensembl version suffix
raw <- raw[setdiff(rownames(raw), "eGFP"), colnames(X)]       # same cells, GFP excluded from size measures

# ---- Ground truth (shared QC: same 745 cells for both readouts) --------------
truth_full <- load_gfp_truth(PROCESSED_DIR, upper_pct = 0.99)        # all cells + qc_pass flag
truth <- truth_full[truth_full$qc_pass & truth_full$Cell_ID %in% colnames(X), ]
meta  <- read.csv(file.path(PROCESSED_DIR, "cell_metadata.csv"), stringsAsFactors = FALSE)
truth$plate <- factor(meta$plate[match(truth$Cell_ID, meta$cell_id)])

# ---- Technical covariates (per cell) ------------------------------------------
# Size factor recovered exactly from raw / normalized counts on shared genes
# (logNormCounts(log = FALSE) divides each cell by its centred size factor).
shared <- intersect(rownames(X), rownames(raw))
covariates <- tibble(
  Cell_ID       = colnames(X),
  size_factor   = Matrix::colSums(raw[shared, ]) / colSums(X[shared, ]),
  total_umis    = Matrix::colSums(raw),
  genes_detected = Matrix::colSums(raw > 0)
)

# ---- Target lists -------------------------------------------------------------
if (!file.exists(TARBASE_RDS))
  stop("tarbase.rds not found at ", TARBASE_RDS,
       "\nRun setup_tools.sh, or set MIRNA_TARBASE_RDS to the file's location.")
tarbase_all <- readRDS(TARBASE_RDS)
# Same filter as miReact's TarBase mode for human: up_down == "DOWN" only.
tarbase_entries <- tarbase_all[tarbase_all$mirna == TARGET_MIRNA &
                               tarbase_all$species == "Homo sapiens" &
                               !is.na(tarbase_all$up_down) & tarbase_all$up_down == "DOWN", ]
# geneId is Ensembl for ~98.5% of entries; this avoids needing miReact's symbol table.
tarbase_entries$gid <- sub("\\.\\d+$", "", as.character(tarbase_entries$geneId))
TARBASE <- intersect(unique(tarbase_entries$gid), rownames(X))

ts_candidates <- file.path(BASE_DIR, c("resources/human/Targets__combined_124-3p_124-3p.2_506-3p.rds",
                                       "resources/Targets__combined_124-3p_124-3p.2_506-3p.rds",
                                       "resources/human/Targets__combined_124-3p_124-3p.2_506-3p_top1000.rds"))
ts_path <- ts_candidates[file.exists(ts_candidates)][1]
if (is.na(ts_path)) stop("No TargetScan target file found under resources/")
TS_FULL_LIST <- !grepl("_top1000", ts_path)
targetscan <- readRDS(ts_path)
targetscan <- targetscan[order(targetscan$Cumulative.weighted.context...score), ]
TS_RANKED  <- unique(targetscan$ensembl_gene_id)                      # best first
TS_TOP200  <- intersect(TS_RANKED[1:N_TOP_TS], rownames(X))
TS_ANY     <- intersect(TS_RANKED, rownames(X))                         # any predicted site

NON_TARGETS <- setdiff(rownames(X), union(TARBASE, TS_ANY))

cat(sprintf("TarBase %s (DOWN): %d genes in matrix (the earlier miReact-based symbol mapping gave 3,017)\n",
            TARGET_MIRNA, length(TARBASE)))
cat(sprintf("TargetScan file: %s%s\n  top %d: %d genes in matrix; any predicted site: %d\n",
            basename(ts_path), if (TS_FULL_LIST) "" else "  (top-1000 only: overlap analyses limited)",
            N_TOP_TS, length(TS_TOP200), length(TS_ANY)))
cat(sprintf("Non-targets (neither list): %d genes\n", length(NON_TARGETS)))

LISTS <- list(`TarBase` = TARBASE, `TargetScan top200` = TS_TOP200)

# ---- Precomputed transforms used by the scoring variants ----------------------
gene_mean <- rowMeans(X)
keep      <- gene_mean > 0
X         <- X[keep, ]; gene_mean <- gene_mean[keep]
LISTS     <- lapply(LISTS, intersect, rownames(X))
NON_TARGETS <- intersect(NON_TARGETS, rownames(X))
LOGRATIO  <- log2((X + PSEUDOCOUNT) / (gene_mean + PSEUDOCOUNT))       # gene-centred log ratio
RANKS     <- apply(LOGRATIO, 2, rank, ties.method = "average")         # within-cell ranks of centred values

# ---- Scoring variants (all: higher = more miRNA activity = targets lower) -----
# A  raw sum                     -- our current method
# B  equal weight per gene       -- each gene divided by its own mean first
# C  mean gene-centred log ratio -- B on log scale (damps outliers)
# C2 C minus non-target mean     -- background-corrected log ratio
# D  mean within-cell rank       -- rank-based, close to miReact's statistic
# (A "D minus non-target ranks" step would be identical to D: within a cell the
#  ranks always sum to the same total, so the non-target mean is fixed by the
#  target mean. So the background step is applied to C instead.)
score_variant <- function(genes, variant) {
  genes <- intersect(genes, rownames(X))
  s <- switch(variant,
    A  = calculate_activity(counts = X, targets = genes),
    B  = { v <- colSums(X[genes, , drop = FALSE] / gene_mean[genes]); -log2(v / mean(v)) + 1 },
    C  = -colMeans(LOGRATIO[genes, , drop = FALSE]),
    C2 = -(colMeans(LOGRATIO[genes, , drop = FALSE]) - colMeans(LOGRATIO[NON_TARGETS, , drop = FALSE])),
    D  = -colMeans(RANKS[genes, , drop = FALSE]),
    stop("Unknown variant ", variant))
  setNames(as.numeric(s), colnames(X))
}
VARIANTS <- c(A = "A: raw sum (current)", B = "B: equal weight per gene", C = "C: mean log ratio",
              C2 = "C2: log ratio minus non-targets", D = "D: mean within-cell rank")

# ---- Evaluation shortcuts -------------------------------------------------------
as_scores <- function(s) tibble(Cell_ID = names(s), score = unname(s))

# Full stratified evaluation with bootstrap CIs (4 rows).
evaluate <- function(s, label, n_boot = N_BOOT)
  evaluate_score_vs_gfp(as_scores(s), "score", label, truth_full, n_boot = n_boot, seed = SEED)

# Fast single rho, no bootstrap (for loops over random draws).
quick_rho <- function(s, readout = "gfp_facs", cell_set = "induced") {
  d <- truth[if (cell_set == "induced") truth$induced else TRUE, ]
  cor(s[d$Cell_ID], d[[readout]], method = "spearman", use = "complete.obs")
}

# Run expr with a fixed seed without disturbing anything else.
seeded <- function(seed, expr) { set.seed(seed); expr }

cat("Setup complete.\n\n")
