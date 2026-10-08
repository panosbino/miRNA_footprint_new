# TarBase diagnostics

Why does our aggregative method do much worse with the TarBase list than ranked
methods (miReact, miTEA-HiRes) given the same list? Each script tests one explanation.

| Script | Question | Main output |
|---|---|---|
| `01_aggregation_ladder.R` | At which step from raw sum to rank-based scoring does TarBase performance recover? | `01_ladder.pdf` |
| `02_weighting.R` | Do a few highly expressed genes dominate the sum? Does removing them help? | `02_concentration.pdf`, `02_drop_top_genes.pdf` |
| `03_list_quality.R` | How many TarBase genes respond? Does evidence type, HEK293 origin or list size matter? | `03_per_gene_rho.pdf`, `03_subsets.pdf`, `03_size_curve.pdf` |
| `04_size_confound.R` | Is the TarBase sum measuring cell size? Does that explain the fluorescence-vs-mRNA gap? | `04_size_confound.pdf`, `04_partial.csv` |

`tarbase_diag_common.R` is shared setup, sourced by every script. Scripts are
independent; run them in any order.

## Running

From the repository root:

```bash
export MIRNA_BASE_DIR="$PWD"
Rscript scripts/HEK_SS3/tarbase_diagnostics/01_aggregation_ladder.R
```

Outputs go to `analysis/HEK_SS3/tarbase_diagnostics/`.

Optional settings (environment variables):

| Variable | Default | Meaning |
|---|---|---|
| `MIRNA_TARBASE_RDS` | `tools/miReact/data/tarbase.rds` | miReact's bundled TarBase file |
| `DIAG_N_BOOT` | 200 | Bootstrap resamples per confidence interval |
| `DIAG_N_RANDOM` | 20 | Random draws for nulls and subset references |

For final figures, use `DIAG_N_BOOT=1000` and `DIAG_N_RANDOM=100`.

## Notes

- **TarBase file.** It ships with miReact, which is not in git. On a fresh clone,
  run `setup_tools.sh` or point `MIRNA_TARBASE_RDS` at a copy.
- **TarBase gene mapping** uses TarBase's own Ensembl `geneId` column, so miReact
  isn't needed. The gene count may differ slightly from the 3,017 mapped in
  `mireact_comparison_rebuild.R` via miReact's symbol table; the setup prints both.
- **TargetScan file.** The full list is used if present, otherwise the committed
  top-1,000 file. With top-1,000 only, "TarBase with a TargetScan site" in script 03
  means "with a top-1,000 site".
- **Ladder step E.** A "mean target rank minus mean non-target rank" step is
  identical to step D: ranks in a cell always sum to the same total. The
  background correction is applied to the log-ratio score (C2) instead.
- **Script 04 residuals.** If cell size is biologically linked to GFP level,
  regressing it out also removes real signal. Read the residual results together
  with the direct covariate correlations.
- Requires the shared helpers in `scripts/Utils.R` from the `gfp-stratified` branch.
