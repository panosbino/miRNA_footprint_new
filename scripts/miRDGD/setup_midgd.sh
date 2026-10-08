#!/usr/bin/env bash
# Set up miDGD (Zamani, Rasmussen et al., bioRxiv 10.64898/2026.05.29.727918, v2) for the
# HEK293 miR-124 comparison. Run from the shared tools folder.
#   bash setup_midgd.sh                 # code, checkpoints, env, model-1 feature lists
#   GET_ZENODO=1 bash setup_midgd.sh    # also fetch the 4.7 GB dataset for models 2-5 feature lists
set -euo pipefail
TOOLS=${TOOLS:-$PWD}
cd "$TOOLS"

# 1. Code + all 7 trained checkpoints, pinned to the commits inspected on 2026-10-07.
#    miDGD_paper holds the checkpoints; miDGD (minimal repo) holds mrna_genes.csv.
#    The inference function learn_new_representation() is identical in both.
[[ -d miDGD_paper ]] || git clone https://github.com/JakobSkouPedersenLab/miDGD_paper.git
git -C miDGD_paper checkout -q 2d249b60ce79a2bfbbbe794503a543e77908f7b1
[[ -d miDGD ]] || git clone https://github.com/JakobSkouPedersenLab/miDGD.git
git -C miDGD checkout -q b7ce3fdbe1dd062c54f3533aa1e08129f1bf7eae

# 2. Inference-only environment. requirement.txt also pulls wandb, scanpy, jupyter: not needed here.
if ! conda env list | grep -q '^midgd '; then
  conda create -y -n midgd python=3.10
  conda run -n midgd pip install "torch>=2.0" --index-url https://download.pytorch.org/whl/cpu
  conda run -n midgd pip install "pandas>=2.0" "numpy>=1.24" "scipy>=1.10" tqdm
fi

# 3. Feature lists for model 1-tcga (18,393 mRNA / 755 miRNA). Order checked against the
#    checkpoint: GMM-mean initialisation gives mean per-miRNA Spearman 0.45 on the bundled
#    TCGA test samples, and collapses to one component with shuffled genes.
mkdir -p midgd_features/1-tcga
cp miDGD/models/mrna_genes.csv midgd_features/1-tcga/mrna.txt
head -1 miDGD/data/tcga_mirna.tsv | tr '\t' '\n' | tail -n +2 > midgd_features/1-tcga/mirna.txt

# 4. Feature lists for models 2-5: column headers of each model's training matrices, read
#    straight from the open Zenodo zip (record 21885644, CC-BY 4.0) without unpacking it.
if [[ "${GET_ZENODO:-0}" == 1 ]]; then
  [[ -f midgd_data.zip ]] || wget -O midgd_data.zip \
    "https://zenodo.org/records/21885644/files/miDGD_paper%20data.zip?download=1"
  echo "dd4e2051bcaad2846316980bd5ee0fd9  midgd_data.zip" | md5sum -c -
  conda run -n midgd python - <<'EOF'
import zipfile, pandas as pd, pathlib, re
# Training files per checkpoint, taken from notebook/*.ipynb in miDGD_paper.
MAP = {
 "1-tcga":                    ("TCGA/miDGD/tcga_mrna.tsv",                        "TCGA/miDGD/tcga_mirna.tsv"),
 "2-gtex":                    ("GTEx/miDGD/gtex_mrna.tsv",                        "GTEx/miDGD/gtex_mirna.tsv"),
 "3a-tcga-gtex":              ("TCGA_GTEx/miDGD/full/tcga_gtex_mrna.tsv",         "TCGA_GTEx/miDGD/full/tcga_gtex_mirna.tsv"),
 "3b-tcga-gtex-collapsed":    ("TCGA_GTEx/miDGD/full/tcga_gtex_mrna_collapsed.tsv","TCGA_GTEx/miDGD/full/tcga_gtex_mirna_collapsed.tsv"),
 "4a-tcga-gtex-r2":           ("TCGA_GTEx_R2/miDGD/TCGA_GTEx_R2_mrna.tsv",        "TCGA_GTEx_R2/miDGD/TCGA_GTEx_R2_mirna.tsv"),
 "4b-tcga-gtex-r2-collapsed": ("TCGA_GTEx_R2/miDGD/TCGA_GTEx_R2_mrna.tsv",        "TCGA_GTEx_R2/miDGD/TCGA_GTEx_R2_mirna_collapsed.tsv"),
 "5-tcga-gtex-r2-smartseq":   ("TCGA_GTEx_SC_R2/miDGD/TCGA_GTEx_SC_R2_mrna.tsv",  "TCGA_GTEx_SC_R2/miDGD/TCGA_GTEx_SC_R2_mirna.tsv"),
}
zf = zipfile.ZipFile("midgd_data.zip")
names = zf.namelist()
def member(rel):
    hits = [n for n in names if n.endswith(rel)]
    if len(hits) != 1: raise SystemExit(f"{rel}: {len(hits)} matches in zip")
    return hits[0]
for model, (mr, mi) in MAP.items():
    out = pathlib.Path("midgd_features") / model
    out.mkdir(parents=True, exist_ok=True)
    for rel, fn in ((mr, "mrna.txt"), (mi, "mirna.txt")):
        cols = pd.read_table(zf.open(member(rel)), index_col=0, nrows=0).columns
        target = out / fn
        if model == "1-tcga" and target.exists():   # cross-check against the GitHub list
            same = target.read_text().split() == list(cols)
            print(f"1-tcga {fn}: GitHub list {'matches' if same else 'DIFFERS FROM'} Zenodo header")
            continue
        target.write_text("\n".join(cols) + "\n")
    mirs = (out / "mirna.txt").read_text().split()
    print(f"{model}: {len((out/'mrna.txt').read_text().split())} mRNA, {len(mirs)} miRNA; "
          f"miR-124 entries: {[m for m in mirs if re.search(r'-mi[Rr]-124(-|$)', m)] or 'none'}")
EOF
fi
echo "Done. Feature lists in $TOOLS/midgd_features/<model>/{mrna,mirna}.txt"
