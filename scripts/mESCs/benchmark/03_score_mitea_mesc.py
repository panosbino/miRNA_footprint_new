"""
mESC benchmark, step 3: miTEA-HiRes with each target set (custom-list route).

Prerequisites
  - scripts/mESCs/export_for_mitea_mesc.R has been run
    (writes analysis/mESCs/mitea_input/counts_data/counts_symbol.txt)
  - step 0 has been run (writes analysis/mESCs/benchmark/target_sets_symbols.tsv)

Run with the mitea environment's own interpreter, inside an allocation:
  /cfs/klemming/projects/supr/naiss2024-6-235/conda/envs/mitea/bin/python3 \
      scripts/mESCs/benchmark/03_score_mitea_mesc.py

Scoring
  Per cell: genes sorted by miTEA's normalised expression (ascending, same call
  as run_mitea_mesc.py), XL-mHG test of target enrichment among the lowest,
  score = -log10(p). The test is called directly with miTEA's own parameters
  (X = mitea_hires.utils._MHG_X_PARAM, L = number of genes). compute_stats_per_cell()
  builds the 0/1 target vector with a list lookup, which is O(targets) per gene
  and takes hours for the 8,000-gene TarBase set; a set lookup gives the same
  vector. A cross-check against compute_stats_per_cell() runs first and stops
  the script if any p-value differs.

  p = 0 is written as NaN (as in run_mitea_mesc.py). The summary step reports
  and drops such cells per (method, target set); it does not hide them.

Output: analysis/mESCs/benchmark/scores_mitea.csv  (method, target_set, cell, score)
"""

import os
from multiprocessing import Pool

import numpy as np
import pandas as pd
import xlmhg
from mitea_hires import process_data
from mitea_hires.utils import compute_stats_per_cell, _MHG_X_PARAM

BASE_DIR = os.environ.get("MIRNA_BASE_DIR", "/cfs/klemming/projects/supr/naiss2024-6-235/miRNA_footprint_new")
DATA_PATH = f"{BASE_DIR}/analysis/mESCs/mitea_input/counts_data/"   # a directory, not a file
BENCH_DIR = f"{BASE_DIR}/analysis/mESCs/benchmark"
N_CPUS = int(os.environ.get("SLURM_CPUS_PER_TASK", "8"))
CHECK_SET = "miRTarBase_high"   # small set, so the slow reference route is quick
N_CHECK_CELLS = 5

counts_norm, _ = process_data(data_path=DATA_PATH, dataset_name="mESC_KO_control_benchmark",
                              data_type="scRNAseq", preprocess=False)
print(f"Loaded: {counts_norm.shape[0]} genes x {counts_norm.shape[1]} cells")

symbols = pd.read_csv(f"{BENCH_DIR}/target_sets_symbols.tsv", sep="\t")
present = set(counts_norm.index)
SETS = {}
for name, grp in symbols.groupby("target_set", sort=False):
    genes = sorted(set(grp["gene_symbol"]) & present)
    SETS[name] = genes
    print(f"{name:22s} {grp.shape[0]:5d} symbols, {len(genes):5d} in miTEA matrix")

# Each cell's ranking is computed once and reused for every set.
RANKED = {cell: list(counts_norm.loc[:, cell].sort_values().index) for cell in counts_norm.columns}


def pval_fast(ranked_list, target_set):
    v = np.uint8([g in target_set for g in ranked_list])
    _, _, pval = xlmhg.xlmhg_test(v, X=_MHG_X_PARAM, L=len(ranked_list))
    return pval


# ---- cross-check against miTEA's own per-cell function ---------------------------------
check_targets = SETS[CHECK_SET]
mti = pd.DataFrame({"miRNA": CHECK_SET, "Target Gene": check_targets})
for cell in list(counts_norm.columns)[:N_CHECK_CELLS]:
    ranked = counts_norm.loc[:, cell].sort_values()
    _, _, ref, _ = compute_stats_per_cell(cell, ranked, [CHECK_SET], mti, debug=False)
    fast = pval_fast(list(ranked.index), set(check_targets))
    if not np.isclose(ref[0], fast, rtol=1e-12, atol=0):
        raise RuntimeError(f"Cross-check failed for {cell}: compute_stats_per_cell p={ref[0]}, direct p={fast}")
print(f"Cross-check passed: direct XL-mHG call reproduces compute_stats_per_cell() on {N_CHECK_CELLS} cells.")


def score_one_set(name):
    target_set = set(SETS[name])
    if not target_set:
        return pd.DataFrame(columns=["method", "target_set", "cell", "score"])
    p = pd.Series({cell: pval_fast(rl, target_set) for cell, rl in RANKED.items()}, dtype=float)
    n_zero = int((p == 0).sum())
    score = (-np.log10(p)).replace([np.inf, -np.inf], np.nan)
    print(f"{name:22s} done: {n_zero} cells with p = 0 (written as NaN)", flush=True)
    return pd.DataFrame({"method": "miTEA-HiRes", "target_set": name, "cell": score.index, "score": score.values})


with Pool(min(N_CPUS, len(SETS))) as pool:
    out = pd.concat(pool.map(score_one_set, list(SETS)), ignore_index=True)

out.to_csv(f"{BENCH_DIR}/scores_mitea.csv", index=False)
print(f"Wrote {BENCH_DIR}/scores_mitea.csv ({out['target_set'].nunique()} target sets)")
