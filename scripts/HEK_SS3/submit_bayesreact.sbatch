#!/bin/bash
#SBATCH --job-name=bayesreact_comparison
#SBATCH --account=naiss2024-6-235
#SBATCH --partition=shared
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --time=04:00:00
#SBATCH --output=bayesreact_comparison_%j.out
#SBATCH --error=bayesreact_comparison_%j.err

# EDIT BEFORE SUBMITTING:
#   --account: confirm this matches your actual allocation (guessed from
#     paths used elsewhere in this project -- not independently verified).
#   --partition: "shared" is a guess -- check `sinfo` or PDC's own docs for
#     the correct partition name/policy on Dardel; this has NOT been
#     confirmed against real Dardel partition names.
#   --time: 4 hours is a rough guess given bayesReact's MCMC fitting is the
#     most compute-intensive step built so far in this project, but this
#     has never actually been timed. If it's not enough, this job will be
#     killed partway through with no partial results saved (bayesReact_core
#     runs as one atomic call) -- consider running interactively first with
#     a small subset to get a real timing estimate before committing to a
#     long unattended job.

module load R  # confirm the exact module name/version on Dardel -- guessed,
                # not verified against this system's actual module list

Rscript /cfs/klemming/projects/snic/naiss2024-6-235/miRNA_footprint_new/scripts/HEK_SS3/comparisons/bayesreact_comparison_rebuild.R
