#!/bin/bash

#SBATCH --job-name=bayesreact_comparison
#SBATCH --account=naiss2026-3-153
#SBATCH --partition=shared
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --time=04:00:00
#SBATCH --output=bayesreact_comparison_%j.out
#SBATCH --error=bayesreact_comparison_%j.err

ml PDC
module load R

Rscript /cfs/klemming/projects/snic/naiss2024-6-235/miRNA_footprint_new/scripts/HEK_SS3/bayesreact_comparison_rebuild.R
