#!/bin/bash
# One array task per unit (row of analysis/PSCSR/benchmark/data/units.tsv).
# Set --array to 1-<number of units> after step 0; --account/--partition/module
# taken from scripts/HEK_SS3/submit_bayesreact.sbatch. Time not measured for
# 2,310 cells -- check the first task before trusting the limit.
#SBATCH --job-name=bayesreact_pscsr
#SBATCH --account=naiss2026-3-153
#SBATCH --partition=shared
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --time=12:00:00
#SBATCH --array=1-30
#SBATCH --output=bayesreact_pscsr_%A_%a.out
#SBATCH --error=bayesreact_pscsr_%A_%a.err

ml PDC
module load R
export LANG=C.UTF-8
export MIRNA_BASE_DIR=/cfs/klemming/projects/supr/naiss2024-6-235/miRNA_footprint_new
Rscript $MIRNA_BASE_DIR/scripts/PSCSR/benchmark/05_run_bayesreact_pscsr.R
