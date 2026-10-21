# ==============================================================================
# PSCSR benchmark, step 6: miReact motif mode (no target list), one run for all
# 7-mers; step 8 reads the rows of the unit motifs. Same call as
# scripts/HEK_SS3/run_mireact_default_launch.R. mireact() submits its own sbatch
# job (run this on a login node); output is already activity-oriented
# (wrapper3.R flips the sign).
#
# miReact writes its own sbatch.script with no account/partition and runs it as
# "Rscript --vanilla", which ignores this session's .libPaths(). Both are fixed
# through environment variables, which sbatch reads and passes into the job:
# SBATCH_ACCOUNT / SBATCH_PARTITION (read by sbatch) and R_LIBS (read by R even
# with --vanilla; needed for Regmex in tools/R_library).
# The job log is slurm-<jobid>.out in miReact's working folder (miReact's message
# says Rscript-<jobid>.out, which is not what SLURM writes).
# ==============================================================================
source(file.path(Sys.getenv("MIRNA_BASE_DIR", "/cfs/klemming/projects/supr/naiss2024-6-235/miRNA_footprint_new"),
                 "scripts/PSCSR/benchmark/config_pscsr.R"))
R_LIB <- file.path(BASE_DIR, "tools/R_library")
.libPaths(c(R_LIB, .libPaths()))
if (!requireNamespace("Regmex", quietly = TRUE)) stop("Regmex is not installed in ", R_LIB, " -- miReact's addMotifs() needs it")
Sys.setenv(R_LIBS = paste(c(R_LIB, Sys.getenv("R_LIBS")), collapse = ":"),
           SBATCH_ACCOUNT = "naiss2026-3-153", SBATCH_PARTITION = "shared")   # as the other submit scripts
source(file.path(MIREACT_DIR, "code/mireact.R"))
for (f in c("hs.seqXmot.utr3_mrs_7mer.rds", "hs.seqXmot.counts.utr3_mrs_7mer.rds"))
  if (!file.exists(file.path(MIREACT_DIR, "motif.models", f))) stop("Missing miReact motif model: ", f)
exp_path <- file.path(DATA_DIR, "scran_normalized_linear.rds")
OUT_FILE <- "mireact_motif_pscsr.rds"
msg <- mireact(exp = exp_path, motifs = 7, species = "hs", seq.type = "utr3", out.file = OUT_FILE,
               out.meonly = TRUE, mail = NULL, install.dir = MIREACT_DIR)
cat(msg, "\n")
wd <- sub("Follow progress in (.*)/Rscript-\\[jobid\\]\\.out", "\\1", msg)
saveRDS(list(wd = wd, out_file = OUT_FILE, submitted_at = Sys.time()), file.path(OUT_DIR, "mireact_motif_job_tracking.rds"))
cat(sprintf("Expected output: %s/%s -- step 8 reads it via mireact_motif_job_tracking.rds\n", wd, OUT_FILE))
cat("Check it was accepted: squeue -u $USER (job name 'Rscript'); log: slurm-<jobid>.out in that folder\n")
