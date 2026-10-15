#!/usr/bin/env Rscript
# ==============================================================================
# miDGD predicted miR-124-3p vs GFP fluorescence, all QC cells and Dox-induced only.
#
# Reads every <model>_mirna_prop.tsv written by run_midgd.py (cells x miRNAs,
# decoder proportions) and scores the miR-124 column with the SHARED
# evaluate_score_vs_gfp() from Utils.R: same QC cells (GFP mRNA < 99th pct),
# same cell sets, same Spearman + percentile-bootstrap CI as every other tool.
#
# Sign: higher predicted miR-124 = more miRNA, so a positive rho is expected,
# matching the "higher = more activity" convention of the other scores.
# ==============================================================================

BASE_DIR      <- Sys.getenv("MIRNA_BASE_DIR", "/cfs/klemming/projects/supr/naiss2024-6-235/miRNA_footprint_new")
PROCESSED_DIR <- file.path(BASE_DIR, "datasets/HEK_SS3/processed")
MIDGD_DIR     <- Sys.getenv("MIDGD_RESULTS_DIR", file.path(BASE_DIR, "analysis/"))
OUT_DIR       <- Sys.getenv("MIRNA_COMPARISON_DIR", file.path(BASE_DIR, "analysis/HEK_SS3/comparisons"))
MIRNA         <- "hsa-miR-124-3p"
GFP_UPPER_PCT <- 0.99
N_BOOT        <- 1000

suppressPackageStartupMessages(library(tidyverse))
source(file.path(BASE_DIR, "scripts/Utils.R"))   # load_gfp_truth(), evaluate_score_vs_gfp()
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

gfp_truth <- load_gfp_truth(PROCESSED_DIR, upper_pct = GFP_UPPER_PCT)
qc_cells  <- gfp_truth$Cell_ID[gfp_truth$qc_pass]

pred_files <- list.files(MIDGD_DIR, pattern = "_mirna_prop\\.tsv$", full.names = TRUE)
if (length(pred_files) == 0) stop("No *_mirna_prop.tsv files in ", MIDGD_DIR)

results <- list()
for (f in pred_files) {
  model <- sub("_mirna_prop\\.tsv$", "", basename(f))
  cat(sprintf("\n== %s ==\n", model))
  pred <- read.delim(f, row.names = 1, check.names = FALSE)

  # miR-124 column: exact name, else a single collapsed 'mir-124' (models 3b/4b/5).
  col <- if (MIRNA %in% colnames(pred)) MIRNA else
    grep("^hsa-mi[Rr]-124(-3p)?$", colnames(pred), value = TRUE)
  if (length(col) != 1) {
    cat(sprintf("  skipped: no unambiguous miR-124 column (found: %s)\n",
                if (length(col)) paste(col, collapse = ", ") else "none"))
    next
  }

  # Cell IDs must be <plate>_<barcode>, as in gfp_correlation_data.csv.
  n_match <- sum(rownames(pred) %in% gfp_truth$Cell_ID)
  cat(sprintf("  %d predicted cells; %d match GFP table; %d of %d QC-passing cells covered\n",
              nrow(pred), n_match, sum(qc_cells %in% rownames(pred)), length(qc_cells)))
  if (n_match == 0) {
    stop("No cell IDs match. Example prediction ID: '", rownames(pred)[1],
         "', example GFP-table ID: '", gfp_truth$Cell_ID[1], "'. ",
         "Rebuild the miDGD input with columns named <plate>_<barcode> ",
         "(bare barcodes can repeat across plates).")
  }
  if (anyDuplicated(rownames(pred))) stop("Duplicate cell IDs in ", basename(f))

  scores <- data.frame(Cell_ID = rownames(pred), midgd = pred[[col]], stringsAsFactors = FALSE)

  # Sanity: miDGD may place HEK293 cells where miR-124 is ~0 and nearly constant.
  s <- scores$midgd[scores$Cell_ID %in% qc_cells]
  cat(sprintf("  %s proportion in QC cells: median %.2e, range [%.2e, %.2e], %d distinct values\n",
              col, median(s), min(s), max(s), length(unique(s))))

  res <- evaluate_score_vs_gfp(scores, "midgd", sprintf("miDGD %s (%s)", model, col),
                               gfp_truth, n_boot = N_BOOT)
  results[[model]] <- res
}

if (length(results) == 0) stop("No model produced a usable miR-124 prediction.")
fluo <- bind_rows(results) |> filter(readout == "fluorescence")
out_file <- file.path(OUT_DIR, "midgd_vs_gfp_fluorescence.csv")
write.csv(fluo, out_file, row.names = FALSE)

cat("\nGFP fluorescence, Spearman rho [95% bootstrap CI]:\n")
print(fluo |> select(method, cell_set, rho, ci_low, ci_high, p_value, n), n = Inf)
cat("Saved:", out_file, "\n")
