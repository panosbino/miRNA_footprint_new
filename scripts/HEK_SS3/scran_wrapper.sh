#!/bin/bash -l
#SBATCH -A naiss2026-3-153
#SBATCH -p shared
#SBATCH -J build_scran_dataset
#SBATCH -t 04:00:00
#SBATCH -n 1
#SBATCH -c 8
#SBATCH -o build_scran_dataset_%j.out
#SBATCH -e build_scran_dataset_%j.err

# Adjust -t (walltime), -c (cores) and add --mem= if your allocation needs an
# explicit memory request -- scran's quickCluster/computeSumFactors step can
# be memory-hungry depending on cell/gene count; these are placeholders.

set -euo pipefail

ml PDC
ml R
ml Rbio

export OPENBLAS_NUM_THREADS=1
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export BLAS_NUM_THREADS=1

Rscript "/cfs/klemming/projects/snic/naiss2024-6-235/miRNA_footprint/scripts/HEK_SS3/build_annotated_scran_dataset.R"
