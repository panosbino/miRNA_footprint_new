"""
Step 3 (all datasets): score every target list with miTEA-HiRes.

    <mitea env>/bin/python3 3_score_mitea.py <mESC|HEK|PSCSR>

Reads <out_dir>/dataset_info.json and target_lists_symbols.tsv (step 1) and the
dataset's miTEA count export (mitea_counts_dir). Per cell, genes are sorted by
miTEA's normalised expression (lowest first) and the XL-mHG test asks whether
the list's genes are enriched among the lowest; score = -log10(p).

The test is called directly with miTEA's own parameters (X = _MHG_X_PARAM,
L = number of genes) and the 0/1 target vector built by array lookup, because
compute_stats_per_cell() looks every gene up in a Python list, which takes hours
for lists of thousands of genes. A cross-check against compute_stats_per_cell()
runs first and stops the script on any difference. p = 0 is written as NaN; the
summary reports such cells. Run by run_benchmark_v3.sbatch (step 'mitea').

Output: <out_dir>/scores/miTEA-HiRes.tsv -- rows list_id, columns cells.
"""
import json
import os
import sys
from multiprocessing import Pool

import numpy as np
import pandas as pd
import xlmhg
from mitea_hires import process_data
from mitea_hires.utils import compute_stats_per_cell, _MHG_X_PARAM

BASE_DIR = os.environ.get("MIRNA_BASE_DIR", "/cfs/klemming/projects/supr/naiss2024-6-235/miRNA_footprint_new")
OUT_DIRS = {"mESC": "analysis/mESCs/benchmark_v3", "HEK": "analysis/HEK_SS3/benchmark_v3", "PSCSR": "analysis/PSCSR/benchmark_v3"}
if len(sys.argv) < 2 or sys.argv[1] not in OUT_DIRS:
    sys.exit("Give the dataset name as the first argument: mESC, HEK or PSCSR")
out_dir = os.path.join(BASE_DIR, OUT_DIRS[sys.argv[1]])
info = json.load(open(os.path.join(out_dir, "dataset_info.json")))
N_CPUS = int(os.environ.get("SLURM_CPUS_PER_TASK", "8"))
N_CHECK_CELLS = 5

counts_norm, _ = process_data(data_path=info["mitea_counts_dir"] + "/", dataset_name=f"benchmark_v3_{info['dataset']}",
                              data_type="scRNAseq", preprocess=False)
print(f"{info['dataset']}: {counts_norm.shape[0]} genes x {counts_norm.shape[1]} cells")

symbols = pd.read_csv(os.path.join(out_dir, "target_lists_symbols.tsv"), sep="\t")
present = set(counts_norm.index)
LISTS = {lid: sorted(set(g["gene_symbol"]) & present) for lid, g in symbols.groupby("list_id")}
LISTS = {k: v for k, v in LISTS.items() if v}
print(f"{len(LISTS)} non-empty lists")

# Each cell's ranking once, reused for every list, stored as gene positions so the
# 0/1 target vector is one array lookup instead of ~10,000 Python lookups. The order
# comes from the same pandas sort as miTEA's own code, so ties stay in the same order.
GENES = list(counts_norm.index)
GENE_POS = {g: i for i, g in enumerate(GENES)}
RANKED = {cell: np.array([GENE_POS[g] for g in counts_norm.loc[:, cell].sort_values().index]) for cell in counts_norm.columns}
TARGET_MASK = {lid: np.isin(np.arange(len(GENES)), [GENE_POS[g] for g in genes]).astype(np.uint8) for lid, genes in LISTS.items()}


def pvalue(ranked_positions, mask):
    is_target = mask[ranked_positions]                     # 0/1 in this cell's gene order
    _, _, p = xlmhg.xlmhg_test(is_target, X=_MHG_X_PARAM, L=len(ranked_positions))
    return p


# Cross-check on the smallest list with >= 10 genes
check_id = min((k for k, v in LISTS.items() if len(v) >= 10), key=lambda k: len(LISTS[k]))
mti = pd.DataFrame({"miRNA": check_id, "Target Gene": LISTS[check_id]})
for cell in list(counts_norm.columns)[:N_CHECK_CELLS]:
    ranked = counts_norm.loc[:, cell].sort_values()
    _, _, ref, _ = compute_stats_per_cell(cell, ranked, [check_id], mti, debug=False)
    fast = pvalue(RANKED[cell], TARGET_MASK[check_id])
    if not np.isclose(ref[0], fast, rtol=1e-12, atol=0):
        raise RuntimeError(f"Cross-check failed for {cell}: compute_stats_per_cell p={ref[0]}, direct p={fast}")
print(f"Cross-check passed on {N_CHECK_CELLS} cells.")


def score_list(list_id):
    mask = TARGET_MASK[list_id]
    p = np.array([pvalue(RANKED[c], mask) for c in counts_norm.columns], dtype=float)
    with np.errstate(divide="ignore"):
        s = -np.log10(p)
    s[~np.isfinite(s)] = np.nan
    return list_id, s


with Pool(min(N_CPUS, len(LISTS))) as pool:
    results = dict(pool.map(score_list, list(LISTS)))

scores = pd.DataFrame.from_dict(results, orient="index", columns=counts_norm.columns)
os.makedirs(os.path.join(out_dir, "scores"), exist_ok=True)
scores.to_csv(os.path.join(out_dir, "scores", "miTEA-HiRes.tsv"), sep="\t", index_label="list_id")
print(f"Wrote scores/miTEA-HiRes.tsv ({scores.shape[0]} lists); cells with p = 0: {int(np.isnan(scores.values).sum())}")
