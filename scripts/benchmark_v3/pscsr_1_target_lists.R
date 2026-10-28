# ==============================================================================
# PSCSR-seq V2, step 1: target lists for every selected miRNA family and its
# controls.
#
# Benchmarked: each family chosen in step 0 (top families per cell line).
# Negative control: miR-122-5p's seed family (not expressed in these lines),
#   rebuilt for EACH benchmarked family, because that family's targets are
#   removed first (CONTROL_EXCLUSION: any database by default).
# Random: N_RANDOM same-size draws per list from genes that are not targets of
#   the benchmarked family (same exclusion).
# Motif methods: each family's site; control = the miR-122 site; random = one
#   shared pool of N_RANDOM random 7-mers (draw k uses motif k for every family).
# Co-expression: correlations pooled WITHIN cell lines; across lines, cell-line
#   differences would dominate every correlation.
#
# Output (analysis/PSCSR/benchmark_v3/): target_lists.rds, target_lists_symbols.tsv, dataset_info.json
# ==============================================================================
source(file.path(Sys.getenv("MIRNA_BASE_DIR", "/cfs/klemming/projects/supr/naiss2024-6-235/miRNA_footprint_new"),
                 "scripts/benchmark_v3/benchmark_config.R"))
source(file.path(BASE_DIR, "scripts/benchmark_v3/benchmark_lib.R"))
ds <- get_dataset("PSCSR")
data_dir <- file.path(ds$out_dir, "data")
set.seed(RANDOM_SEED)

# ---- Data from step 0 ---------------------------------------------------------------------------
expr  <- readRDS(ds$expression)
cells <- read_tsv(file.path(data_dir, "cells.tsv"), show_col_types = FALSE)
fams  <- read_tsv(file.path(data_dir, "families.tsv"), show_col_types = FALSE) |> mutate(members = strsplit(members, ";"))
stopifnot(identical(colnames(expr), cells$cell))
utr <- load_utr_table("hs")
universe <- intersect(rownames(expr)[rowMeans(expr) > 0], utr$gid)
cat(sprintf("Scored genes (with a 3'UTR): %d; benchmarked families: %d\n", length(universe), nrow(fams)))

control <- human_family_of(ds$control_mirna, read_family_info())
ctrl_name <- sub("^hsa-", "", ds$control_mirna); ctrl_site <- revcomp_dna(control$seed)
if (ctrl_site %in% fams$site) stop("The control miRNA shares a seed with a benchmarked family")

# ---- Evidence for all families at once --------------------------------------------------------------
evidence <- load_human_evidence(families = c(setNames(fams$members, fams$family), setNames(list(control$members), ctrl_name)),
                                seeds = c(setNames(fams$seed, fams$family), setNames(control$seed, ctrl_name)),
                                universe = universe, utr = utr)

detected <- rownames(expr)[rowMeans(expr == 0) < COEXP_MAX_ZERO_FRAC]
line_groups <- split(cells$cell, cells$cell_line)
coexp <- function(genes) coexpression_filter(genes, expr, line_groups, detected)

# ---- Lists per family ------------------------------------------------------------------------------------
rows <- list()
for (f in fams$family) {
  cat(sprintf("Building lists for %s ...\n", f))
  bench_lists <- build_target_lists(evidence[[f]], coexp, with_mirdip = TRUE)
  excluded <- benchmark_exclusion(evidence[[f]])
  ctrl_lists  <- build_target_lists(remove_targets(evidence[[ctrl_name]], excluded), coexp, with_mirdip = TRUE)
  rows[[f]]   <- assemble_list_rows(f, bench_lists, ctrl_lists, build_random_lists(bench_lists, excluded, universe))
}
lists <- bind_rows(rows)

# ---- Motif lists ------------------------------------------------------------------------------------------------
random_pool <- random_motifs(N_RANDOM, avoid = c(fams$site, ctrl_site))
motif_lists <- bind_rows(
  tibble(family = fams$family, kind = "benchmark", draw = NA_integer_, motifs = as.list(fams$site)),
  tibble(family = fams$family, kind = "negative_control", draw = NA_integer_, motifs = list(ctrl_site)),
  expand_grid(family = fams$family, draw = seq_len(N_RANDOM)) |> mutate(kind = "random", motifs = as.list(random_pool[draw])))

families <- bind_rows(
  tibble(family = fams$family, role = "benchmark", name = fams$family, sites = as.list(fams$site), members = fams$members),
  tibble(family = NA_character_, role = "control", name = ctrl_name, sites = list(ctrl_site), members = list(control$members)))
save_target_lists(ds, lists, universe, families, motif_lists, utr, extra = list(random_motif_pool = random_pool))
cat(sprintf("\nWrote %s/target_lists.rds\n", ds$out_dir))
