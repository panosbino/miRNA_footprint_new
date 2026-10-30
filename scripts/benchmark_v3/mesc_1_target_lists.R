# ==============================================================================
# mESC, step 1: target lists for the benchmarked families and their controls.
#
# Benchmarked (all lost in the Drosha KO):
#   "miR-17/291/292 (combined)"  the three seed families together, as in v1/v2
#   "miR-17", "miR-291a", "miR-292a"  each family alone, on its full target lists
#       (targets shared with the other two families are kept)
# Negative control: miR-122-5p (not expressed in mESCs), same rules, after removing
#   the targets of ALL THREE families -- the KO removes all three, so a control that
#   kept e.g. miR-17 targets would separate KO from Control for the miR-291 test.
# Random: N_RANDOM same-size draws from genes that are not targets of any of the
#   three families (the dataset's control_exclusion decides what counts as a target).
# Motif methods: combined = the three sites summed; single family = its site;
#   control = the miR-122 site; random = motifs from a pool of N_RANDOM random 7-mers
#   (three per draw for the combined set, one for a single family).
#
# Output (analysis/mESCs/benchmark_v3/): target_lists.rds, target_lists_symbols.tsv, dataset_info.json
# ==============================================================================
source(file.path(Sys.getenv("MIRNA_BASE_DIR", "/cfs/klemming/projects/supr/naiss2024-6-235/miRNA_footprint_new"),
                 "scripts/benchmark_v3/benchmark_config.R"))
source(file.path(BASE_DIR, "scripts/benchmark_v3/benchmark_lib.R"))
ds <- get_dataset("mESC")
set.seed(RANDOM_SEED)

# ---- Expression and gene universe ---------------------------------------------------
pheno <- read.csv(file.path(BASE_DIR, "datasets/mESCs/phenotype.csv"), fileEncoding = "UTF-8-BOM")
pheno$Exp[pheno$Exp == "WT"] <- "Control"
expr <- read_count_matrix(ds$expression)[, pheno$cells]
utr <- load_utr_table("mm")
universe <- intersect(rownames(expr)[rowMeans(expr) > 0], utr$gid)
cat(sprintf("Expression: %d genes x %d cells; scored genes (expressed, with a 3'UTR): %d\n",
            nrow(expr), ncol(expr), length(universe)))

# ---- Families and seed sites -------------------------------------------------------------
families <- ds$benchmark_families
combined <- names(families)[1]                       # the three families together
singles  <- names(families)[-1]
ctrl     <- ds$control_family
mature   <- read_fasta(MOUSE_DB$mature_fasta)
sites    <- map(families, \(f) unname(vapply(f$seed_mirnas, \(m) site_of(mature[[m]]), "")))
ctrl_site <- site_of(mature[[ctrl$seed_mirnas]])
stopifnot(identical(sites[[combined]], c("GCACTTT", "AGCACTT", "GGCACTT")))   # sites of the original mESC benchmark
for (f in names(sites)) cat(sprintf("%-28s site(s) %s\n", f, paste(sites[[f]], collapse = "+")))
cat(sprintf("%-28s site %s (negative control)\n", ctrl$name, ctrl_site))

# ---- Evidence ----------------------------------------------------------------------------------
evidence <- map(families, \(f) load_mouse_evidence(f, universe, utr))
ctrl_evidence <- load_mouse_evidence(ctrl, universe, utr)
empty_dbs <- names(ctrl_evidence)[map_int(ctrl_evidence, nrow) == 0 & names(ctrl_evidence) != "mirDIP"]
if (length(empty_dbs)) cat(sprintf("*** No %s evidence in: %s -- rerun scripts/extract_targets_mesc.R (it must list %s)\n",
                                   ctrl$name, paste(empty_dbs, collapse = ", "), paste(ctrl$table_labels, collapse = ", ")))

# ---- Co-expression rule: correlations among Control cells --------------------------------------------
detected <- rownames(expr)[rowMeans(expr == 0) < COEXP_MAX_ZERO_FRAC]
control_cells <- pheno$cells[pheno$Exp == "Control"]
coexp <- function(genes) coexpression_filter(genes, expr, list(control_cells), detected)

# ---- Controls: built once, excluding the targets of all three families -------------------------------------
excluded   <- benchmark_exclusion(evidence[[combined]], ds$control_exclusion)
ctrl_lists <- build_target_lists(remove_targets(ctrl_evidence, excluded), coexp, with_mirdip = FALSE)

# ---- Benchmarked lists per family ------------------------------------------------------------------------------
rows <- list()
for (f in names(families)) {
  cat(sprintf("Building lists for %s ...\n", f))
  bench_lists <- build_target_lists(evidence[[f]], coexp, with_mirdip = FALSE)
  rows[[f]] <- assemble_list_rows(f, bench_lists, ctrl_lists, build_random_lists(bench_lists, excluded, universe))
}
lists <- bind_rows(rows)

# ---- Motif lists ---------------------------------------------------------------------------------------------------
random_pool <- random_motifs(N_RANDOM, avoid = c(unlist(sites), ctrl_site))
motif_lists <- map_dfr(names(families), function(f) {
  k <- length(sites[[f]])
  bind_rows(
    tibble(family = f, kind = "benchmark", draw = NA_integer_, motifs = list(sites[[f]])),
    tibble(family = f, kind = "negative_control", draw = NA_integer_, motifs = list(ctrl_site)),
    tibble(family = f, kind = "random", draw = seq_len(N_RANDOM),
           motifs = if (k == 1) as.list(random_pool) else map(seq_len(N_RANDOM), \(i) sample(random_pool, k))))
})

families_tbl <- bind_rows(
  tibble(family = names(families), role = "benchmark", name = names(families), sites = unname(sites)),
  tibble(family = NA_character_, role = "control", name = ctrl$name, sites = list(ctrl_site)))
save_target_lists(ds, lists, universe, families_tbl, motif_lists, utr, extra = list(random_motif_pool = random_pool))
cat(sprintf("\nWrote %s/target_lists.rds\n", ds$out_dir))
