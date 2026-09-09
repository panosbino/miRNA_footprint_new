"""
Run miTEA-HiRes on scran-normalized data, restricted to hsa-miR-124-3p,
both native miRTarBase and TarBase-matched. Same logic as the original
build -- conda env (mitea, Python 3.8) confirmed intact, no changes needed
there. Only paths updated for the new directory structure.

Activate the conda environment first:
    conda activate mitea
"""

import pandas as pd
import numpy as np
from mitea_hires import process_data, compute_mir_activity
from mitea_hires.utils import compute_stats_per_cell

BASE_DIR = "/cfs/klemming/projects/supr/naiss2024-6-235/miRNA_footprint_new"
OUT_DIR = f"{BASE_DIR}/analysis/HEK_SS3/comparisons"
MITEA_INPUT_DIR = f"{OUT_DIR}/mitea_input"

DATA_PATH = f"{MITEA_INPUT_DIR}/counts_data/"   # directory, not file -- process_data() globs this
RESULTS_PATH = f"{MITEA_INPUT_DIR}/results/"
TARGET_MIRNA = "hsa-miR-124-3p"

counts_norm, counts_raw = process_data(
    data_path=DATA_PATH,
    dataset_name="miR124_HEK293_scran_comparison",
    data_type="scRNAseq",
    preprocess=False,   # avoid the 10,000-cell auto-subsample threshold, same as before
)
print(f"Loaded: {counts_norm.shape[0]} genes x {counts_norm.shape[1]} cells")


def score_and_check(pvals, label):
    n_zero = (pvals == 0).sum()
    print(f"*** SANITY CHECK ({label}): {n_zero} of {len(pvals)} cells have p=0 for {TARGET_MIRNA} ***")
    if n_zero == len(pvals):
        raise RuntimeError(f"[{label}] ALL cells have p=0 -- no targets found. Do not proceed.")
    elif n_zero > 0:
        print(f"WARNING ({label}): {n_zero} cells have p=0 -- set to NaN below, not inf.")
    activity_score = -np.log10(pvals.astype(float))
    return activity_score.replace([np.inf, -np.inf], np.nan)


# --- (1) Native miRTarBase scoring -------------------------------------------
# cpus explicitly capped to the REAL interactive allocation (8), not left
# as the default None. Very likely cause of the earlier BrokenPipeError /
# ForkPoolWorker-90 crash: cpus=None probably falls back to something like
# os.cpu_count(), which on a shared HPC node reports the PHYSICAL node's
# total core count, not what's actually allocated to this session --
# causing massive oversubscription (many more workers spawned than real
# cores available) and processes getting killed under resource pressure.
N_CPUS = 8   # MATCH YOUR ACTUAL ALLOCATION -- if you request more/fewer
             # cores for a future run, update this to match, don't leave
             # it out of sync with the real session/job allocation.

miR_list, miR_activity_pvals = compute_mir_activity(
    counts_norm, results_path=RESULTS_PATH, miR_list=[TARGET_MIRNA],
    species="homo_sapiens", debug=True, cpus=N_CPUS,
)
mirtarbase_score = score_and_check(miR_activity_pvals.loc[TARGET_MIRNA], "miRTarBase")
out_mirtarbase = pd.DataFrame({"Cell_ID": mirtarbase_score.index, "mitea_activity": mirtarbase_score.values})
out_mirtarbase.to_csv(f"{MITEA_INPUT_DIR}/mitea_activity_scores_mirtarbase_scran.csv", index=False)
print(f"Wrote mitea_activity_scores_mirtarbase_scran.csv ({out_mirtarbase['mitea_activity'].notna().sum()} non-NaN)")

# --- (2) TarBase scoring, bypassing built-in target loading -----------------
with open(f"{MITEA_INPUT_DIR}/tarbase_124_targets.txt") as f:
    tarbase_targets = [line.strip() for line in f if line.strip()]

tarbase_targets_in_data = [g for g in tarbase_targets if g in counts_norm.index]
print(f"TarBase targets: {len(tarbase_targets)} total, {len(tarbase_targets_in_data)} present in expression matrix")

mti_data_tarbase = pd.DataFrame({
    "miRNA": [TARGET_MIRNA] * len(tarbase_targets_in_data),
    "Target Gene": tarbase_targets_in_data,
})

tarbase_pvals = {}
for cell in counts_norm.columns:
    ranked = counts_norm.loc[:, cell].sort_values()
    _, _, pvals_row, _ = compute_stats_per_cell(cell, ranked, [TARGET_MIRNA], mti_data_tarbase, debug=False)
    tarbase_pvals[cell] = pvals_row[0]

tarbase_pvals = pd.Series(tarbase_pvals)
tarbase_score = score_and_check(tarbase_pvals, "TarBase")
out_tarbase = pd.DataFrame({"Cell_ID": tarbase_score.index, "mitea_activity_tarbase": tarbase_score.values})
out_tarbase.to_csv(f"{MITEA_INPUT_DIR}/mitea_activity_scores_tarbase_scran.csv", index=False)
print(f"Wrote mitea_activity_scores_tarbase_scran.csv ({out_tarbase['mitea_activity_tarbase'].notna().sum()} non-NaN)")
