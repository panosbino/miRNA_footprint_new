"""Predict per-cell miRNA expression from mRNA counts with a pretrained miDGD checkpoint.

Uses the authors' own learn_new_representation() unchanged (GMM-mean initialisation, then
50 AdamW steps on z with mRNA NB loss + GMM loss; decoder and GMM frozen).

Input  --counts : RAW exon UMI counts, genes x cells, TSV, Ensembl gene IDs as row names,
                  ALL genes (not the scran-filtered 10,785-gene matrix).
Output <out>_mirna_prop.tsv : cells x miRNAs, decoder proportions (rows sum to 1)
       <out>_latent.tsv     : cells x latent dims

Example:
  conda run -n midgd python run_midgd.py --repo tools/miDGD_paper \
      --model tools/miDGD_paper/models/4a-tcga-gtex-r2.pth \
      --features tools/midgd_features/4a-tcga-gtex-r2 \
      --counts hek_ss3_exon_umi.tsv --out results/midgd_4a
"""
import argparse, sys, numpy as np, pandas as pd, torch

ap = argparse.ArgumentParser()
ap.add_argument("--repo", required=True, help="path to miDGD_paper clone (provides base/)")
ap.add_argument("--model", required=True)
ap.add_argument("--features", required=True, help="dir with mrna.txt and mirna.txt for this model")
ap.add_argument("--counts", required=True)
ap.add_argument("--out", required=True)
ap.add_argument("--mirna", default="hsa-miR-124-3p")
ap.add_argument("--epochs", type=int, default=50)   # authors' default
ap.add_argument("--seed", type=int, default=1312)
a = ap.parse_args()

sys.path.insert(0, a.repo)  # checkpoints pickle classes from base.*
from torch.utils.data import DataLoader
from base.data.combined import GeneExpressionDatasetCombined
from base.engine.predict import learn_new_representation, device

torch.manual_seed(a.seed)
dgd = torch.load(a.model, map_location=device, weights_only=False)
dgd.eval()

genes = open(f"{a.features}/mrna.txt").read().split()
mirnas = open(f"{a.features}/mirna.txt").read().split()
dec = dgd.decoder
if (len(genes), len(mirnas)) != (dec.n_out_features_mrna, dec.n_out_features_mirna):
    sys.exit(f"Feature lists ({len(genes)}, {len(mirnas)}) do not match checkpoint "
             f"({dec.n_out_features_mrna}, {dec.n_out_features_mirna}): wrong --features dir.")
if a.mirna not in mirnas:
    sys.exit(f"{a.mirna} is not an output of this model. 124 entries: "
             f"{[m for m in mirnas if '124' in m]}")

counts = pd.read_table(a.counts, index_col=0)
counts.index = counts.index.str.replace(r"\.\d+$", "", regex=True)   # drop Ensembl version
if counts.index.duplicated().any():
    counts = counts.groupby(level=0).sum()
if (counts.values % 1 != 0).any():
    sys.exit("Counts are not integers: miDGD needs raw counts (NB likelihood), not normalised values.")

present = counts.index.intersection(genes)
frac = len(present) / len(genes)
print(f"{len(present)} / {len(genes)} model genes present ({100*frac:.1f}%)")
if frac < 0.9:
    print("WARNING: library size is summed over present genes only, but the decoder's softmax "
          "spans all model genes, so many absent genes mis-scale the NB mean.")

# Absent genes -> NaN: skipped by the NB likelihood (same reindex the authors use for
# Smart-seq-total and PSCSR-seq). A gene that is annotated but unobserved stays 0.
X = counts.T.reindex(columns=genes).astype(float)
mirna_dummy = pd.DataFrame(np.nan, index=X.index, columns=mirnas)  # miRNA unobserved; its loss is not optimised
anno = pd.DataFrame({c: "query" for c in ["primary_site", "tissue_type", "cancer_type", "batch"]},
                    index=X.index).assign(color="#000000")

ds = GeneExpressionDatasetCombined(X, mirna_dummy, anno, scaling_type="sum")
loader = DataLoader(ds, batch_size=256, shuffle=False)
rep = learn_new_representation(dgd, loader, test_epochs=a.epochs, learning_rates=1e-2)

with torch.no_grad():
    z = rep.z.detach()
    prop_mirna, _ = dgd.decoder(z.to(device))

# Proportions, deliberately NOT multiplied by a library size: the authors scale by the
# observed miRNA library, which does not exist here, and an mRNA-derived scale would
# re-introduce sequencing depth (correlated with GFP in this dataset).
pd.DataFrame(prop_mirna.cpu().numpy(), index=X.index, columns=mirnas).to_csv(f"{a.out}_mirna_prop.tsv", sep="\t")
pd.DataFrame(z.cpu().numpy(), index=X.index).to_csv(f"{a.out}_latent.tsv", sep="\t")
print(f"Wrote {a.out}_mirna_prop.tsv ({a.mirna} is one column) and {a.out}_latent.tsv")
