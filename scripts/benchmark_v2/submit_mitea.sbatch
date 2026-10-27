#!/bin/bash
# miTEA-HiRes for all three datasets, one array task per dataset, run in parallel:
#   sbatch scripts/benchmark_v2/submit_mitea.sbatch
# Only some datasets: --array=2 (PSCSR only) or --array=0,1 (mESC and HEK).
# Needs step 1 for each dataset, and its miTEA count export (mESC: export_for_mitea_mesc.R;
# HEK: export_for_mitea_rebuild.R; PSCSR: pscsr_0_prepare_data.R).
#SBATCH --job-name=mitea_v2
#SBATCH --account=naiss2026-3-153
#SBATCH --partition=shared
#SBATCH --array=0-2
#SBATCH --nodes=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --time=24:00:00
#SBATCH --output=mitea_v2_%A_%a.out
#SBATCH --error=mitea_v2_%A_%a.err

DATASETS=(mESC HEK PSCSR)
DS=${DATASETS[$SLURM_ARRAY_TASK_ID]}
export MIRNA_BASE_DIR=/cfs/klemming/projects/supr/naiss2024-6-235/miRNA_footprint_new
export LANG=C.UTF-8
PY=/cfs/klemming/projects/supr/naiss2024-6-235/conda/envs/mitea/bin/python3

echo "miTEA-HiRes on $DS with $SLURM_CPUS_PER_TASK cores, started $(date)"
$PY $MIRNA_BASE_DIR/scripts/benchmark_v2/3_score_mitea.py "$DS"
echo "finished $(date)"
