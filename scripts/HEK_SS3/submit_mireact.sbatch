#!/bin/bash
#SBATCH --job-name=mireact_comparison
#SBATCH --account=naiss2026-3-153
#SBATCH --partition=shared
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --time=01:00:00
#SBATCH --output=mireact_comparison_%j.out
#SBATCH --error=mireact_comparison_%j.err

# SAME UNVERIFIED-GUESS CAVEATS AS submit_bayesreact.sbatch -- confirm
# --account, --partition, and the module name below against your actual
# Dardel setup before submitting; none of these have been checked.
#
# --time=01:00:00 is a much shorter estimate than bayesReact's, since this
# is wcmod.p() (a closed-form rank-sum statistic) rather than MCMC fitting
# -- genuinely lighter computation, not just an arbitrary smaller guess.
# Still never actually timed, though. --cpus-per-task=8 must match
# N_WORKERS in the R script itself (currently also 8) -- if you change one,
# change the other, they don't sync automatically.
ml PDC
module load R 

Rscript /cfs/klemming/projects/snic/naiss2024-6-235/miRNA_footprint_new/scripts/HEK_SS3/mireact_comparison_rebuild.R
