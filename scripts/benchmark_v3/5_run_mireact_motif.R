# ==============================================================================
# Step 5 (all datasets): miReact motif mode, run inside the current job.
#
#   Rscript 5_run_mireact_motif.R <mESC|HEK|PSCSR>
#
# One run scores every 7-mer, so the benchmarked, control and random motifs are
# all read from the same output. mireact() would submit its own SLURM job; here
# the working folder is prepared the same way (runparameters.Rdata and e.Rdata in
# a folder inside the miReact install directory) and miReact's wrapper3.R is run
# directly, so it works inside the pipeline job (run_benchmark_v3.sbatch).
# R_LIBS points to tools/R_library, which holds Regmex.
# Output is activity-oriented (wrapper3.R flips the sign).
# Writes <out_dir>/mireact_motif_job.rds (working folder and output file), which
# the summary step reads.
# ==============================================================================
source(file.path(Sys.getenv("MIRNA_BASE_DIR", "/cfs/klemming/projects/supr/naiss2024-6-235/miRNA_footprint_new"),
                 "scripts/benchmark_v3/benchmark_config.R"))
ds <- dataset_from_args()
.libPaths(c(R_LIBRARY, .libPaths()))
if (!requireNamespace("Regmex", quietly = TRUE)) stop("Regmex is not installed in ", R_LIBRARY)
Sys.setenv(R_LIBS = paste(c(R_LIBRARY, Sys.getenv("R_LIBS")), collapse = ":"))

tl <- readRDS(file.path(ds$out_dir, "target_lists.rds"))
out_file <- sprintf("mireact_motif_%s.rds", ds$name)

# Working folder inside the install directory (wrapper3.R is called as ../code/wrapper3.R)
wd <- file.path(MIREACT_DIR, sprintf("o_benchmark_v3_%s_%s", ds$name, format(Sys.time(), "%Y%m%d_%H%M%S")))
dir.create(wd)
runparameters <- list(exp = tl$expression, motifs = 7, species = DATASETS[[ds$name]]$species, seq.type = "utr3",
                      out.file = out_file, out.meonly = TRUE, verbose = FALSE, install.dir = MIREACT_DIR,
                      mail = NULL, memory = "128g", wt = "48:00:00", tarbaserun = FALSE, wd = wd)
save(runparameters, file = file.path(wd, "runparameters.Rdata"))
file.symlink(tl$expression, file.path(wd, "e.Rdata"))            # as mireact() does for a file path

cat(sprintf("miReact motif mode for %s in %s (started %s)\n", ds$name, wd, format(Sys.time())))
owd <- setwd(wd)                                                  # wrapper3.R expects to run inside the folder
status <- system2("Rscript", c("--vanilla", "../code/wrapper3.R"),
                  stdout = "wrapper3.log", stderr = "wrapper3.log")   # R_LIBS is inherited from this session
setwd(owd)

result <- file.path(wd, out_file)
if (status != 0 || !file.exists(result))
  stop(sprintf("miReact motif mode failed (exit status %s). See %s", status, file.path(wd, "wrapper3.log")))
saveRDS(list(wd = wd, out_file = out_file, finished = Sys.time()), file.path(ds$out_dir, "mireact_motif_job.rds"))
cat(sprintf("Done (%s). Result: %s\n", format(Sys.time()), result))
