# ==============================================================================
# mESC, step 1: target lists for the benchmarked families and their controls.
#
# Benchmarked: miR-17, miR-291 and miR-292 seed families together (the miRNAs
#   lost in the Drosha KO that dominate mESCs), as in the original benchmark.
# Negative control: miR-122-5p targets, built with the same rules after removing
#   the benchmarked families' targets (CONTROL_EXCLUSION: any database by default).
# Random: N_RANDOM same-size draws from genes that are not targets of the
#   benchmarked families (same exclusion).
# Motif methods (bayesReact, miReact motif mode): the three seed sites summed;
#   control = the miR-122 site; random = 3 motifs drawn from a pool of N_RANDOM
#   random 7-mers.
#
# Output (analysis/mESCs/benchmark_v2/): target_lists.rds,
#   target_lists_symbols.tsv, dataset_info.json
# ==============================================================================
source(file.path(Sys.getenv("MIRNA_BASE_DIR", "/cfs/klemming/projects/supr/naiss2024-6-235/miRNA_footprint_new"),
                 "scripts/benchmark_v2/benchmark_config.R"))
source(file.path(BASE_DIR, "scripts/benchmark_v2/benchmark_lib.R"))
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

# ---- Families, seed sites ---------------------------------------------------------------
bench_name <- names(ds$benchmark_families)
bench <- ds$benchmark_families[[1]]
ctrl  <- ds$control_family
mature <- read_fasta(MOUSE_DB$mature_fasta)
bench_sites <- vapply(bench$seed_mirnas, \(m) site_of(mature[[m]]), "")
ctrl_site   <- site_of(mature[[ctrl$seed_mirnas]])
stopifnot(identical(unname(bench_sites), c("GCACTTT", "AGCACTT", "GGCACTT")))   # sites used in the original mESC benchmark
cat(sprintf("Seed sites: %s = %s; control %s = %s\n", bench_name, paste(bench_sites, collapse = "+"), ctrl$name, ctrl_site))

# ---- Evidence from the databases ------------------------------------------------------------
bench_evidence <- load_mouse_evidence(bench, universe, utr)
ctrl_evidence  <- load_mouse_evidence(ctrl, universe, utr)
empty_dbs <- names(ctrl_evidence)[map_int(ctrl_evidence, nrow) == 0 & names(ctrl_evidence) != "mirDIP"]
if (length(empty_dbs)) cat(sprintf("*** No %s evidence for %s in: %s -- rerun scripts/extract_targets_mesc.R (it must list %s)\n",
                                   ctrl$name, ctrl$name, paste(empty_dbs, collapse = ", "), paste(ctrl$table_labels, collapse = ", ")))

# ---- Co-expression rule: correlations among Control cells --------------------------------------
# KO cells lack mature miRNAs; with them included, the KO-vs-Control shift itself
# would create the correlations the benchmark then scores.
detected <- rownames(expr)[rowMeans(expr == 0) < COEXP_MAX_ZERO_FRAC]
control_cells <- pheno$cells[pheno$Exp == "Control"]
coexp <- function(genes) coexpression_filter(genes, expr, list(control_cells), detected)

# ---- Target lists ---------------------------------------------------------------------------------------
bench_lists <- build_target_lists(bench_evidence, coexp, with_mirdip = FALSE)
excluded <- benchmark_exclusion(bench_evidence)
ctrl_lists  <- build_target_lists(remove_targets(ctrl_evidence, excluded), coexp, with_mirdip = FALSE)
random_lists <- build_random_lists(bench_lists, excluded, universe)
lists <- assemble_list_rows(bench_name, bench_lists, ctrl_lists, random_lists)

# ---- Motif lists for motif-based methods -----------------------------------------------------------------
random_pool <- random_motifs(N_RANDOM, avoid = c(bench_sites, ctrl_site))
motif_lists <- bind_rows(
  tibble(family = bench_name, kind = "benchmark", draw = NA_integer_, motifs = list(unname(bench_sites))),
  tibble(family = bench_name, kind = "negative_control", draw = NA_integer_, motifs = list(ctrl_site)),
  tibble(family = bench_name, kind = "random", draw = seq_len(N_RANDOM),
         motifs = map(seq_len(N_RANDOM), \(i) sample(random_pool, length(bench_sites)))))

families <- tibble(family = bench_name, role = c("benchmark", "control"),
                   name = c(bench_name, ctrl$name), sites = list(unname(bench_sites), ctrl_site))
save_target_lists(ds, lists, universe, families, motif_lists, utr,
                  extra = list(random_motif_pool = random_pool))
cat(sprintf("\nWrote %s/target_lists.rds\n", ds$out_dir))
