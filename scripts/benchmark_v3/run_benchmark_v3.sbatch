#!/bin/bash
# ==============================================================================
# Benchmark v3: the whole pipeline for all three datasets in one job.
#
#   sbatch scripts/benchmark_v3/run_benchmark_v3.sbatch
#
# Per dataset (mESC, then HEK, then PSCSR):
#   prepare        PSCSR only: matrices and the top families per line (pscsr_0)
#   lists          target lists with negative and random controls (<dataset>_1)
#   then IN PARALLEL, each with its own share of the cores:
#     score          our 5 variants + miReact on the lists (2)
#     mitea          miTEA-HiRes (3)
#     bayesreact     bayesReact, every motif (4)
#     mireact_motif  miReact motif mode (5)
#     midgd          miDGD, HEK and PSCSR only (6 + run_midgd.py)
#   summary        tables, heatmaps and per-family plots (7)
#
# Choose datasets / steps with environment variables, e.g. to redo only the
# summaries:  sbatch --export=ALL,STEPS="summary" scripts/benchmark_v3/run_benchmark_v3.sbatch
#             sbatch --export=ALL,DATASETS="PSCSR",STEPS="mitea summary" ...
# Each step writes logs/benchmark_v3/<dataset>_<step>.log; this job's own log
# lists every step as done / FAILED, and a failed step does not stop the others
# (only a failed "lists" skips the rest of that dataset).
#
# Before the first run: the human database files (benchmark_config.R, HUMAN_DB),
# scripts/extract_targets_mesc.R run once with miR-122, and the miTEA count exports
# for mESC and HEK (export_for_mitea_mesc.R, export_for_mitea_rebuild.R).
# Check the time limit allowed on your partition; steps can be resumed with STEPS.
# ==============================================================================
#SBATCH --job-name=benchmark_v3
#SBATCH --account=naiss2026-3-153
#SBATCH --partition=shared
#SBATCH --nodes=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=180G
#SBATCH --time=72:00:00
#SBATCH --output=benchmark_v3_%j.out
#SBATCH --error=benchmark_v3_%j.err

# ---- Settings --------------------------------------------------------------------------------
export MIRNA_BASE_DIR=${MIRNA_BASE_DIR:-/cfs/klemming/projects/supr/naiss2024-6-235/miRNA_footprint_new}
DATASETS=${DATASETS:-"mESC HEK PSCSR"}
STEPS=${STEPS:-"prepare lists score mitea bayesreact mireact_motif midgd summary"}
MITEA_PY=${MITEA_PY:-/cfs/klemming/projects/supr/naiss2024-6-235/conda/envs/mitea/bin/python3}
MIDGD_PY=${MIDGD_PY:-/cfs/klemming/projects/supr/naiss2024-6-235/conda/envs/midgd/bin/python}   # check this path
# cores per parallel step (sum <= --cpus-per-task)
CORES_SCORE=8; CORES_MITEA=10; CORES_BAYES=8; CORES_MIDGD=4          # miReact motif mode uses 1

ml PDC
module load R
export LANG=C.UTF-8
export R_LIBS=$MIRNA_BASE_DIR/tools/R_library${R_LIBS:+:$R_LIBS}
cd "$MIRNA_BASE_DIR"
S=scripts/benchmark_v3
LOG_DIR=logs/benchmark_v3/job_${SLURM_JOB_ID:-local}
mkdir -p "$LOG_DIR"

out_dir() { case $1 in mESC) echo analysis/mESCs/benchmark_v3 ;; HEK) echo analysis/HEK_SS3/benchmark_v3 ;; PSCSR) echo analysis/PSCSR/benchmark_v3 ;; esac; }
lists_script() { case $1 in mESC) echo mesc_1_target_lists.R ;; HEK) echo hek_1_target_lists.R ;; PSCSR) echo pscsr_1_target_lists.R ;; esac; }
want() { [[ " $STEPS " == *" $1 "* ]]; }
stamp() { date "+%Y-%m-%d %H:%M:%S"; }

# run <dataset> <step> <command...>: log to its own file, record done / FAILED
run() {
  local ds=$1 step=$2; shift 2
  local log="$LOG_DIR/${ds}_${step}.log"
  echo "[$(stamp)] $ds $step: started (log: $log)"
  if "$@" > "$log" 2>&1; then
    echo "done" > "$LOG_DIR/${ds}_${step}.status"; echo "[$(stamp)] $ds $step: done"
  else
    echo "FAILED" > "$LOG_DIR/${ds}_${step}.status"; echo "[$(stamp)] $ds $step: FAILED -- see $log"
  fi
}

midgd_step() {
  local ds=$1 out counts
  out=$(out_dir "$ds")
  if [ ! -x "$MIDGD_PY" ]; then echo "miDGD python not found: $MIDGD_PY (set MIDGD_PY)"; return 1; fi
  if [ "$ds" = HEK ]; then
    Rscript $S/6_export_midgd_counts_hek.R || return 1
    counts=$out/midgd/hek_exon_umi.tsv
  else
    counts=$out/data/midgd_counts.tsv
  fi
  mkdir -p "$out/midgd"
  OMP_NUM_THREADS=$CORES_MIDGD "$MIDGD_PY" scripts/miRDGD/run_midgd.py --repo tools/miDGD_paper \
    --model tools/miDGD_paper/models/4a-tcga-gtex-r2.pth --features tools/midgd_features/4a-tcga-gtex-r2 \
    --counts "$counts" --out "$out/midgd/4a"
}

echo "[$(stamp)] benchmark v3: datasets [$DATASETS], steps [$STEPS], $SLURM_CPUS_PER_TASK cores; logs in $LOG_DIR"

for ds in $DATASETS; do
  echo "==================== $ds ===================="
  if want prepare && [ "$ds" = PSCSR ]; then run "$ds" prepare Rscript $S/pscsr_0_prepare_data.R; fi
  if want lists; then
    run "$ds" lists Rscript $S/$(lists_script "$ds")
    if [ "$(cat "$LOG_DIR/${ds}_lists.status")" = FAILED ]; then echo "Skipping the rest of $ds."; continue; fi
  fi

  # independent steps, in parallel
  if want score;         then run "$ds" score env N_CORES=$CORES_SCORE Rscript $S/2_score_target_lists.R "$ds" & fi
  if want mitea;         then run "$ds" mitea env SLURM_CPUS_PER_TASK=$CORES_MITEA "$MITEA_PY" $S/3_score_mitea.py "$ds" & fi
  if want bayesreact;    then run "$ds" bayesreact env SLURM_CPUS_PER_TASK=$CORES_BAYES Rscript $S/4_run_bayesreact.R "$ds" & fi
  if want mireact_motif; then run "$ds" mireact_motif Rscript $S/5_run_mireact_motif.R "$ds" & fi
  if want midgd && [ "$ds" != mESC ]; then run "$ds" midgd midgd_step "$ds" & fi
  wait

  if want summary; then run "$ds" summary Rscript $S/7_summarise.R "$ds"; fi
done

echo "==================== overview ===================="
for f in "$LOG_DIR"/*.status; do [ -e "$f" ] && printf "%-28s %s\n" "$(basename "$f" .status)" "$(cat "$f")"; done
if grep -l FAILED "$LOG_DIR"/*.status > /dev/null 2>&1; then echo "Some steps FAILED -- rerun them with STEPS=... after fixing."; fi
echo "[$(stamp)] finished"
