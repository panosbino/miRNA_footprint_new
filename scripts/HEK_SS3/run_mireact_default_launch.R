library(tidyverse)

# ---------------------------------------------------------------------------
# Launches miReact's DEFAULT motif-search mode -- finally usable now that
# the required precomputed files (hs.seqXmot.utr3_mrs_7mer.rds and the
# .counts. variant) are confirmed present in tools/miReact/motif.models/.
#
# RUN THIS SCRIPT DIRECTLY (Rscript, on the login node) -- do NOT wrap it in
# your own sbatch submission. mireact() submits its own SLURM job
# internally (confirmed by reading mireact.R directly) and returns almost
# immediately; the actual computation happens asynchronously afterward.
#
# IMPORTANT: mireact() creates a TIMESTAMPED working directory
# (o[currentDateAndTime]) that can't be predicted in advance. Its return
# value tells us exactly where -- this script captures and saves that path
# so collect_mireact_default_results.R (run LATER, once the job finishes)
# knows where to look, rather than guessing/searching for the newest
# timestamped folder.
# ---------------------------------------------------------------------------

BASE_DIR <- "/cfs/klemming/projects/snic/naiss2024-6-235/miRNA_footprint_new"
TOOLS_DIR <- file.path(BASE_DIR, "tools")
MIREACT_DIR <- file.path(TOOLS_DIR, "miReact")
PROCESSED_DIR <- file.path(BASE_DIR, "datasets/HEK_SS3/processed")
OUT_DIR <- file.path(BASE_DIR, "analysis/HEK_SS3/comparisons")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

.libPaths(c(file.path(TOOLS_DIR, "R_library"), .libPaths()))

source(file.path(MIREACT_DIR, "code/mireact.R"))

# Confirm the required files are actually where mireact() expects them,
# BEFORE submitting anything -- fail loudly here, not deep inside an async
# SLURM job that's much more annoying to debug.
required_files <- c("hs.seqXmot.utr3_mrs_7mer.rds", "hs.seqXmot.counts.utr3_mrs_7mer.rds")
for (f in required_files) {
  full_path <- file.path(MIREACT_DIR, "motif.models", f)
  if (!file.exists(full_path)) stop(sprintf("Missing required file: %s", full_path))
  cat(sprintf("Confirmed present: %s (%.1f MB)\n", f, file.info(full_path)$size / 1e6))
}

# --- Point directly at the existing scran-normalized RDS file --------------
# Passed as a PATH (not loaded into memory) -- mireact() handles this via
# try(load(...)) with a readRDS() fallback (verified directly in
# wrapper3.R), so either an in-memory matrix or a file path works
# correctly. Using the path avoids loading the whole matrix into this
# lightweight launcher session at all.
exp_path <- file.path(PROCESSED_DIR, "scran_normalized_linear.rds")
stopifnot(file.exists(exp_path))

OUT_FILE <- "mireact_default_hek293_scran.rds"

result_message <- mireact(
  exp = exp_path,
  motifs = 7,          # default, stated explicitly for clarity
  species = "hs",       # default, stated explicitly for clarity
  seq.type = "utr3",    # default, stated explicitly for clarity
  out.file = OUT_FILE,
  out.meonly = TRUE,    # default -- motif-activity matrix only
  mail = NULL,           # unset -- no mail server assumed available
  install.dir = MIREACT_DIR
)

cat("\n", result_message, "\n")

# Parse the working directory out of the returned message and save it, so
# the results-collection script (run later) knows exactly where to look.
wd_path <- sub("Follow progress in (.*)/Rscript-\\[jobid\\]\\.out", "\\1", result_message)
tracking_info <- list(wd = wd_path, out_file = OUT_FILE, submitted_at = Sys.time())
saveRDS(tracking_info, file.path(OUT_DIR, "mireact_default_job_tracking.rds"))

cat(sprintf("\nJob directory: %s\n", wd_path))
cat(sprintf("Expected output file: %s/%s\n", wd_path, OUT_FILE))
cat("Tracking info saved to:", file.path(OUT_DIR, "mireact_default_job_tracking.rds"), "\n")
cat("\nCheck progress with: squeue -u $USER\n")
cat(sprintf("Or watch the log directly: tail -f %s/Rscript-*.out\n", wd_path))
cat("\nOnce the job completes, run collect_mireact_default_results.R (NOT this script again).\n")
