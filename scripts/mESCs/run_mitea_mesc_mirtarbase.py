"""
miTEA-HiRes on the mESC Control/KO data, using its NATIVE bundled mouse
miRTarBase (mmu_MTI_filtered.csv, loaded via load_mti('mus_musculus')),
per miRNA family and combined.

Prerequisite: export_for_mitea_mesc.R has already been run (it writes
mitea_input/counts_data/counts_symbol.txt).
Run with the mitea environment's own interpreter, e.g.
  /cfs/klemming/projects/supr/naiss2024-6-235/conda/envs/mitea/bin/python3 run_mitea_mesc_mirtarbase.py

DESIGN DECISIONS
1. miRTarBase lists targets per INDIVIDUAL miRNA (e.g. mmu-miR-17-5p), not per
   family. A family's target set here = the UNION of the targets of its member
   miRNAs; "combined" = the union across all three families. One XL-mHG test
   per (cell, set), the same structure used for the other tools.
2. Targets are matched to the expression matrix CASE-INSENSITIVELY. The
   bundled file mixes capitalization (a few entries are all-uppercase, human
   style) while mouse symbols are title case, and miTEA-HiRes matches names
   case-sensitively -- a plain match would silently drop those genes.
3. Scoring uses compute_stats_per_cell() directly (same code the native
   pipeline calls per cell), because the native wrapper takes single miRNAs
   and cannot take a family-level union. A cross-check below confirms the
   wiring reproduces the native compute_mir_activity() result exactly for
   one miRNA.

COVERAGE WARNING (measured from the bundled file, not assumed): validated
targets for these miRNAs are few -- roughly 48 unique genes for MIR17, 11 for
MIR291 and 2 for MIR292. XL-mHG is nearly powerless with a handful of genes,
so MIR291 and especially MIR292 results should be read as unreliable. Target
counts are written to mirtarbase_target_counts.csv for the R step to report.
"""

import os
import numpy as np
import pandas as pd
from mitea_hires import process_data, compute_mir_activity
from mitea_hires.utils import compute_stats_per_cell, load_mti

BASE_DIR = "/cfs/klemming/projects/supr/naiss2024-6-235/miRNA_footprint_new"
OUT_DIR = f"{BASE_DIR}/analysis/mESCs"
MITEA_INPUT_DIR = f"{OUT_DIR}/mitea_input"
DATA_PATH = f"{MITEA_INPUT_DIR}/counts_data/"      # a directory, not a file
RESULTS_PATH = f"{MITEA_INPUT_DIR}/results/"
os.makedirs(RESULTS_PATH, exist_ok=True)

N_CPUS = 8                 # MATCH YOUR ALLOCATION (cpus=None over-subscribed the node before)
MIN_TARGETS_WARN = 10
SPECIES = "mus_musculus"   # confirmed in the package's constants

# Family membership: my mapping from the TargetScan family file names to
# mouse miRNA names -- REVIEW THIS. Other names present in the database
# (e.g. mmu-miR-290a-3p, 467b/f/g, the -5p arms of 291a/292a) are deliberately
# NOT included because they are not in these TargetScan families / have
# different seeds.
FAMILIES = {
    "MIR17":  ["mmu-miR-17-5p", "mmu-miR-20a-5p", "mmu-miR-20b-5p",
               "mmu-miR-93-5p", "mmu-miR-106a-5p", "mmu-miR-106b-5p"],
    "MIR291": ["mmu-miR-291a-3p", "mmu-miR-291b-3p", "mmu-miR-294-3p",
               "mmu-miR-295-3p", "mmu-miR-302a-3p", "mmu-miR-302b-3p",
               "mmu-miR-302c-3p", "mmu-miR-302d-3p"],
    "MIR292": ["mmu-miR-292a-3p", "mmu-miR-467a-5p"],
}
FAMILIES["combined"] = sorted({m for v in list(FAMILIES.values()) for m in v})

counts_norm, counts_raw = process_data(
    data_path=DATA_PATH, dataset_name="mESC_KO_control_mirtarbase",
    data_type="scRNAseq", preprocess=False,
)
print(f"Loaded: {counts_norm.shape[0]} genes x {counts_norm.shape[1]} cells")

mti_all = load_mti(SPECIES)
print(f"Native mouse miRTarBase: {mti_all.shape[0]} rows, {mti_all['miRNA'].nunique()} miRNAs")

# case-insensitive lookup from the expression matrix's gene names
idx_upper = {}
for g in counts_norm.index:
    idx_upper.setdefault(str(g).upper(), []).append(g)


def targets_for(mirnas):
    present = [m for m in mirnas if m in set(mti_all["miRNA"])]
    missing = [m for m in mirnas if m not in present]
    sub = mti_all[mti_all["miRNA"].isin(present)]
    raw = sorted(set(sub["Target Gene"].astype(str)))
    matched = sorted({x for t in raw for x in idx_upper.get(t.upper(), [])})
    return present, missing, len(sub), len(raw), matched


def score_targets(matched, label):
    mti_custom = pd.DataFrame({"miRNA": [label] * len(matched), "Target Gene": matched})
    pvals = {}
    for cell in counts_norm.columns:
        ranked = counts_norm.loc[:, cell].sort_values()
        _, _, pv, _ = compute_stats_per_cell(cell, ranked, [label], mti_custom, debug=False)
        pvals[cell] = pv[0]
    return pd.Series(pvals).astype(float)


summary_rows = []
for fam, mirnas in FAMILIES.items():
    present, missing, n_rows, n_raw, matched = targets_for(mirnas)
    print(f"\n=== {fam} ===")
    print(f"member miRNAs found in database: {len(present)}/{len(mirnas)}"
          + (f"  (missing: {', '.join(missing)})" if missing else ""))
    print(f"MTI rows: {n_rows}; unique target genes: {n_raw}; present in expression matrix: {len(matched)}")
    summary_rows.append({"variant": fam, "n_member_mirnas_found": len(present),
                         "n_mti_rows": n_rows, "n_unique_targets": n_raw,
                         "n_targets_in_data": len(matched)})
    if len(matched) == 0:
        print(f"*** SKIPPED {fam}: no targets present in the expression matrix ***")
        continue
    if len(matched) < MIN_TARGETS_WARN:
        print(f"*** WARNING: only {len(matched)} targets (< {MIN_TARGETS_WARN}) -- XL-mHG has almost no "
              f"power here; treat this result as unreliable. ***")

    pvals = score_targets(matched, fam)
    if (pvals == 0).all():
        raise RuntimeError(f"[{fam}] ALL cells have p=0 -- no targets found. Do not proceed.")
    score = (-np.log10(pvals)).replace([np.inf, -np.inf], np.nan)
    n_distinct = score.nunique()
    print(f"-log10(p): min={score.min():.3f} median={score.median():.3f} max={score.max():.3f}; "
          f"distinct values={n_distinct} of {len(score)} cells")
    if n_distinct < 0.5 * len(score):
        print(f"*** NOTE: only {n_distinct} distinct score values across {len(score)} cells -- the score is "
              f"heavily quantized (expected with very few targets), which limits how well it can separate groups. ***")
    pd.DataFrame({"Cell_ID": score.index, "activity": score.values}).to_csv(
        f"{MITEA_INPUT_DIR}/mitea_scores_mirtarbase_{fam}.csv", index=False)

pd.DataFrame(summary_rows).to_csv(f"{MITEA_INPUT_DIR}/mirtarbase_target_counts.csv", index=False)
print(f"\nWrote per-variant score CSVs and mirtarbase_target_counts.csv to {MITEA_INPUT_DIR}")

# --- Cross-check: bypass wiring vs the NATIVE compute_mir_activity ----------
# One miRNA, exact-case targets, both routes. Expected: identical p-values.
# Wrapped so a failure here (e.g. multiprocessing trouble) cannot block the
# results above.
CHECK_MIR = "mmu-miR-17-5p"
try:
    _, native_p = compute_mir_activity(counts_norm, results_path=RESULTS_PATH,
                                       miR_list=[CHECK_MIR], species=SPECIES, cpus=N_CPUS)
    native = native_p.loc[CHECK_MIR].astype(float)
    exact_targets = sorted(set(mti_all.loc[mti_all["miRNA"] == CHECK_MIR, "Target Gene"].astype(str)))
    exact_in_data = [g for g in exact_targets if g in set(counts_norm.index)]
    bypass = score_targets(exact_in_data, CHECK_MIR)
    common = native.index.intersection(bypass.index)
    max_diff = float(np.nanmax(np.abs(native.loc[common].values - bypass.loc[common].values)))
    print(f"\n*** CROSS-CHECK vs native compute_mir_activity ({CHECK_MIR}, exact-case targets): "
          f"max |p difference| = {max_diff:.3e} over {len(common)} cells ***")
    if max_diff > 1e-9:
        print("*** WARNING: bypass does NOT reproduce the native result -- do not trust the results above until resolved. ***")
    n_fold = len(targets_for([CHECK_MIR])[4])
    print(f"(case-insensitive matching found {n_fold} targets vs {len(exact_in_data)} exact-case for {CHECK_MIR})")
except Exception as e:
    print(f"\n*** CROSS-CHECK could not run: {type(e).__name__}: {e} ***\n"
          f"Main results above were still written; rerun the check separately (try a smaller N_CPUS).")
