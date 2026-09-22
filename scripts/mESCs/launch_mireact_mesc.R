library(tidyverse)

# ---------------------------------------------------------------------------
# Launches miReact's default motif-search mode for mESC data. Same async
# SLURM pattern as the HEK293 default-mode work -- run this directly
# (login node), NOT wrapped in your own sbatch job; mireact() submits its
# own job internally and returns almost immediately.
#
# Only ONE miReact run needed here, not three -- unlike bayesReact,
# miReact's motif-search computes ALL 16,384 7-mers in one pass, so all
# three family motifs (GCACTTT/AGCACTT/GGCACTT) can be extracted as rows
# from a single output matrix. Confirmed mouse motif-model files already
# present in tools/miReact/motif.models/.
# ---------------------------------------------------------------------------

BASE_DIR <- "/cfs/klemming/projects/supr/naiss2024-6-235/miRNA_footprint_new"
TOOLS_DIR <- file.path(BASE_DIR, "tools")
MIREACT_DIR <- file.path(TOOLS_DIR, "miReact")
OUT_DIR <- file.path(BASE_DIR, "analysis/mESCs")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

.libPaths(c(file.path(TOOLS_DIR, "R_library"), .libPaths()))
source(file.path(MIREACT_DIR, "code/mireact.R"))

required_files <- c("mm.seqXmot.utr3_mrs_7mer.rds", "mm.seqXmot.counts.utr3_mrs_7mer.rds")
for (f in required_files) {
  full_path <- file.path(MIREACT_DIR, "motif.models", f)
  if (!file.exists(full_path)) stop(sprintf("Missing required file: %s", full_path))
  cat(sprintf("Confirmed present: %s (%.1f MB)\n", f, file.info(full_path)$size / 1e6))
}

setwd(BASE_DIR)
pheno <- read.delim("./datasets/mESCs/phenotype.csv", sep = ",")
counts <- readRDS("datasets/mESCs/final_counts.rds")
counts <- counts[, colnames(counts) %in% pheno$cells]

# Save prepped expression data for mireact() to read as a file path (same
# try(load())-then-readRDS() fallback confirmed in wrapper3.R -- either an
# in-memory matrix or a saved RDS path works correctly).
exp_path <- file.path(OUT_DIR, "mesc_exp_for_mireact.rds")
saveRDS(as.matrix(counts), exp_path)

OUT_FILE <- "mireact_default_mesc.rds"

result_message <- mireact(
  exp = exp_path,
  motifs = 7,
  species = "mm",
  seq.type = "utr3",
  out.file = OUT_FILE,
  out.meonly = TRUE,
  mail = NULL,
  install.dir = MIREACT_DIR
)

cat("\n", result_message, "\n")

wd_path <- sub("Follow progress in (.*)/Rscript-\\[jobid\\]\\.out", "\\1", result_message)
tracking_info <- list(wd = wd_path, out_file = OUT_FILE, submitted_at = Sys.time())
saveRDS(tracking_info, file.path(OUT_DIR, "mireact_mesc_job_tracking.rds"))

cat(sprintf("\nJob directory: %s\n", wd_path))
cat(sprintf("Expected output file: %s/%s\n", wd_path, OUT_FILE))
cat("\nCheck progress with: squeue -u $USER\n")
cat(sprintf("Or watch the log directly: tail -f %s/Rscript-*.out\n", wd_path))
cat("\nOnce the job completes, run collect_mireact_mesc_results.R (NOT this script again).\n")
