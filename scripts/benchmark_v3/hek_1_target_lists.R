# ==============================================================================
# HEK293, step 1: target lists for the miR-124-3p seed family and its controls.
#
# Benchmarked: all human miRNAs sharing miR-124-3p's seed (nt 2-8), per
#   TargetScan's miR_Family_Info.
# Negative control: miR-122-5p's seed family (not expressed in HEK293), same
#   rules, after removing miR-124's targets (the dataset's control_exclusion in benchmark_config.R).
# Random: N_RANDOM same-size draws from genes that are not miR-124 targets (same exclusion).
# Motif methods: the miR-124 site; control = the miR-122 site; random = N_RANDOM
#   random 7-mers. miDGD: the family's predicted proportion; control = miR-122's
#   (random miRNAs are drawn from miDGD's model in the summary step).
# Target databases: the human bulk files in HUMAN_DB (shared with PSCSR).
#
# Output (analysis/HEK_SS3/benchmark_v3/): target_lists.rds,
#   target_lists_symbols.tsv, dataset_info.json
# ==============================================================================
source(file.path(Sys.getenv("MIRNA_BASE_DIR", "/cfs/klemming/projects/supr/naiss2024-6-235/miRNA_footprint_new"),
                 "scripts/benchmark_v3/benchmark_config.R"))
source(file.path(BASE_DIR, "scripts/benchmark_v3/benchmark_lib.R"))
ds <- get_dataset("HEK")
set.seed(RANDOM_SEED)

# ---- Expression and gene universe -------------------------------------------------------
expr <- read_count_matrix(ds$expression)
gfp  <- read.csv(file.path(BASE_DIR, "datasets/HEK_SS3/processed/gfp_correlation_data.csv"))
utr  <- load_utr_table("hs")
universe <- intersect(rownames(expr)[rowMeans(expr) > 0], utr$gid)
cat(sprintf("Expression: %d genes x %d cells; scored genes (with a 3'UTR): %d\n", nrow(expr), ncol(expr), length(universe)))

# ---- Families -----------------------------------------------------------------------------------
family_info <- read_family_info()
bench <- human_family_of(ds$benchmark_mirnas, family_info)
ctrl  <- human_family_of(ds$control_mirna, family_info)
bench_name <- sub("^hsa-", "", ds$benchmark_mirnas); ctrl_name <- sub("^hsa-", "", ds$control_mirna)
bench_site <- revcomp_dna(bench$seed); ctrl_site <- revcomp_dna(ctrl$seed)
stopifnot(bench_site == "GTGCCTT")                       # the miR-124 site used in every earlier HEK analysis
cat(sprintf("%s family: %s (site %s)\n%s family: %s (site %s)\n", bench_name, paste(bench$members, collapse = ", "), bench_site,
            ctrl_name, paste(ctrl$members, collapse = ", "), ctrl_site))

# ---- Evidence ----------------------------------------------------------------------------------------
evidence <- load_human_evidence(families = setNames(list(bench$members, ctrl$members), c(bench_name, ctrl_name)),
                                seeds = setNames(c(bench$seed, ctrl$seed), c(bench_name, ctrl_name)),
                                universe = universe, utr = utr)

# ---- Co-expression rule: correlations among Dox-induced cells -------------------------------------------
# With uninduced cells included, the induced-vs-uninduced contrast (which drives
# GFP, the ground truth) would create the correlations.
detected <- rownames(expr)[rowMeans(expr == 0) < COEXP_MAX_ZERO_FRAC]
induced_cells <- gfp$cell_id[gfp$Dox_Group == "Dox Induced"]
coexp <- function(genes) coexpression_filter(genes, expr, list(induced_cells), detected)

# ---- Target lists -------------------------------------------------------------------------------------------
bench_lists  <- build_target_lists(evidence[[bench_name]], coexp, with_mirdip = TRUE)
excluded <- benchmark_exclusion(evidence[[bench_name]], ds$control_exclusion)
ctrl_lists   <- build_target_lists(remove_targets(evidence[[ctrl_name]], excluded), coexp, with_mirdip = TRUE)
random_lists <- build_random_lists(bench_lists, excluded, universe)
lists <- assemble_list_rows(bench_name, bench_lists, ctrl_lists, random_lists)

# ---- Motif lists ---------------------------------------------------------------------------------------------
random_pool <- random_motifs(N_RANDOM, avoid = c(bench_site, ctrl_site))
motif_lists <- bind_rows(
  tibble(family = bench_name, kind = "benchmark", draw = NA_integer_, motifs = list(bench_site)),
  tibble(family = bench_name, kind = "negative_control", draw = NA_integer_, motifs = list(ctrl_site)),
  tibble(family = bench_name, kind = "random", draw = seq_len(N_RANDOM), motifs = as.list(random_pool)))

families <- tibble(family = bench_name, role = c("benchmark", "control"), name = c(bench_name, ctrl_name),
                   sites = list(bench_site, ctrl_site), members = list(bench$members, ctrl$members))
save_target_lists(ds, lists, universe, families, motif_lists, utr, extra = list(random_motif_pool = random_pool))
cat(sprintf("\nWrote %s/target_lists.rds\n", ds$out_dir))
