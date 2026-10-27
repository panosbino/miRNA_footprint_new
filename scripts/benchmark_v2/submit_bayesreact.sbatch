#!/bin/bash
# One task per motif:  sbatch --array=1-<n> submit_bayesreact.sbatch <mESC|HEK|PSCSR>
# (n is printed by: Rscript scripts/benchmark_v2/4_run_bayesreact.R <dataset>; Ctrl-C after the first line)
# Runtime per motif was not measured for 2,310 cells (PSCSR) -- check the first task.
#SBATCH --job-name=bayesreact_v2
#SBATCH --account=naiss2026-3-153
#SBATCH --partition=shared
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --time=12:00:00
#SBATCH --output=bayesreact_v2_%A_%a.out
#SBATCH --error=bayesreact_v2_%A_%a.err
ml PDC
module load R
export LANG=C.UTF-8
export MIRNA_BASE_DIR=/cfs/klemming/projects/supr/naiss2024-6-235/miRNA_footprint_new
Rscript $MIRNA_BASE_DIR/scripts/benchmark_v2/4_run_bayesreact.R "$1"
