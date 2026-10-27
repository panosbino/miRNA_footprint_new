# ==============================================================================
# Step 5 (all datasets): miReact motif mode -- one run scores every 7-mer, so
# the benchmarked, control and random motifs are all read from the same output.
#
#   Rscript 5_launch_mireact_motif.R <mESC|HEK|PSCSR>     (on a login node)
#
# mireact() writes and submits its own job. Its job script has no account or
# partition and runs "Rscript --vanilla", which ignores this session's library
# path; both are passed as environment variables, which sbatch reads and hands
# to the job (SBATCH_ACCOUNT, SBATCH_PARTITION, R_LIBS -- Regmex is needed).
# Output is activity-oriented (wrapper3.R flips the sign). The job's log is
# slurm-<jobid>.out in miReact's working folder.
# Writes <out_dir>/mireact_motif_job.rds (working folder and output file name).
# ==============================================================================
source(file.path(Sys.getenv("MIRNA_BASE_DIR", "/cfs/klemming/projects/supr/naiss2024-6-235/miRNA_footprint_new"),
                 "scripts/benchmark_v2/benchmark_config.R"))
ds <- dataset_from_args()
.libPaths(c(R_LIBRARY, .libPaths()))
if (!requireNamespace("Regmex", quietly = TRUE)) stop("Regmex is not installed in ", R_LIBRARY)
Sys.setenv(R_LIBS = paste(c(R_LIBRARY, Sys.getenv("R_LIBS")), collapse = ":"),
           SBATCH_ACCOUNT = SLURM_ACCOUNT, SBATCH_PARTITION = SLURM_PARTITION)
source(file.path(MIREACT_DIR, "code/mireact.R"))
tl <- readRDS(file.path(ds$out_dir, "target_lists.rds"))
out_file <- sprintf("mireact_motif_%s.rds", ds$name)
msg <- mireact(exp = tl$expression, motifs = 7, species = DATASETS[[ds$name]]$species, seq.type = "utr3",
               out.file = out_file, out.meonly = TRUE, mail = NULL, install.dir = MIREACT_DIR)
wd <- sub("Follow progress in (.*)/Rscript-\\[jobid\\]\\.out", "\\1", msg)
saveRDS(list(wd = wd, out_file = out_file, submitted = Sys.time()), file.path(ds$out_dir, "mireact_motif_job.rds"))
cat(sprintf("Submitted. Working folder: %s\nCheck: squeue -u $USER (job 'Rscript'); log: slurm-<jobid>.out in that folder.\nResult: %s/%s\n", wd, wd, out_file))
