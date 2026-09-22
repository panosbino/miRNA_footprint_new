"""
Run miTEA-HiRes on mESC data, combined and negative-control target lists,
via the custom bypass (same mechanism as the human TarBase/TargetScan
work). No native mouse miRTarBase needed.

Activate the conda environment first: conda activate mitea
"""

import pandas as pd
import numpy as np
from mitea_hires import process_data
from mitea_hires.utils import compute_stats_per_cell

BASE_DIR = "/cfs/klemming/projects/supr/naiss2024-6-235/miRNA_footprint_new"
OUT_DIR = f"{BASE_DIR}/analysis/mESCs"
MITEA_INPUT_DIR = f"{OUT_DIR}/mitea_input"

DATA_PATH = f"{MITEA_INPUT_DIR}/counts_data/"  # directory, not file
N_CPUS = 8  # match your actual allocation

counts_norm, counts_raw = process_data(
    data_path=DATA_PATH,
    dataset_name="mESC_KO_control",
    data_type="scRNAseq",
    preprocess=False,
)
print(f"Loaded: {counts_norm.shape[0]} genes x {counts_norm.shape[1]} cells")


def score_custom_targets(target_file, label, mirna_label="mESC_target"):
    with open(f"{MITEA_INPUT_DIR}/{target_file}") as f:
        targets = [line.strip() for line in f if line.strip()]
    targets_in_data = [g for g in targets if g in counts_norm.index]
    print(f"{label}: {len(targets)} total, {len(targets_in_data)} present in expression matrix")

    mti_data_custom = pd.DataFrame({
        "miRNA": [mirna_label] * len(targets_in_data),
        "Target Gene": targets_in_data,
    })

    pvals = {}
    for cell in counts_norm.columns:
        ranked = counts_norm.loc[:, cell].sort_values()
        _, _, pvals_row, _ = compute_stats_per_cell(cell, ranked, [mirna_label], mti_data_custom, debug=False)
        pvals[cell] = pvals_row[0]
    pvals = pd.Series(pvals)

    n_zero = (pvals == 0).sum()
    print(f"*** SANITY CHECK ({label}): {n_zero} of {len(pvals)} cells have p=0 ***")
    if n_zero == len(pvals):
        raise RuntimeError(f"[{label}] ALL cells have p=0 -- no targets found. Do not proceed.")

    activity_score = -np.log10(pvals.astype(float))
    activity_score = activity_score.replace([np.inf, -np.inf], np.nan)
    out_df = pd.DataFrame({"Cell_ID": activity_score.index, "activity": activity_score.values})
    out_path = f"{MITEA_INPUT_DIR}/mitea_scores_{label}.csv"
    out_df.to_csv(out_path, index=False)
    print(f"Wrote {out_path} ({out_df['activity'].notna().sum()} non-NaN)")


score_custom_targets("combined_targets.txt", "combined")
score_custom_targets("negctrl_targets.txt", "negctrl")
