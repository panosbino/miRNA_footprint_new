# ==============================================================================
# Why are TarBase / miRTarBase lists empty? (human datasets: HEK, PSCSR)
#
#   Rscript diagnose_target_databases.R <HEK|PSCSR>
#
# Repeats, stage by stage, what load_human_evidence() does for these two
# databases and prints how many rows survive each stage, so the stage where
# they drop to zero is visible. Reads target_lists.rds from step 1 (for the
# miRNA families and the scored genes); changes nothing.
# ==============================================================================
source(file.path(Sys.getenv("MIRNA_BASE_DIR", "/cfs/klemming/projects/supr/naiss2024-6-235/miRNA_footprint_new"),
                 "scripts/benchmark_v2/benchmark_config.R"))
source(file.path(BASE_DIR, "scripts/benchmark_v2/benchmark_lib.R"))
ds <- dataset_from_args()
if (DATASETS[[ds$name]]$species != "hs") stop("This check is for the human datasets (HEK, PSCSR)")

tl <- readRDS(file.path(ds$out_dir, "target_lists.rds"))
universe <- tl$universe
utr <- load_utr_table("hs")
members <- unique(unlist(tl$families$members))
cat(sprintf("%s: %d scored genes; %d family members, e.g. %s\n", ds$name, length(universe), length(members),
            paste(head(members, 5), collapse = ", ")))

cat("\n== Genes per list in target_lists.rds (benchmark / negative control) ==\n")
print(tl$lists |> filter(kind != "random", grepl("^(TarBase|miRTarBase)", list_type)) |>
        group_by(list_type, kind) |> summarise(families = n(), median_genes = median(n_genes), max_genes = max(n_genes), .groups = "drop") |>
        pivot_wider(names_from = kind, values_from = c(median_genes, max_genes)), n = Inf)

stage <- function(label, n) cat(sprintf("  %-55s %10s\n", label, format(n, big.mark = ",")))
similar_names <- function(found) {
  core <- unique(sub("^hsa-(miR|let)-([0-9a-z]+).*$", "\\2", members))[1:3]
  hits <- unique(grep(paste0("-(", paste(core, collapse = "|"), ")([^0-9]|$)"), found, value = TRUE))
  cat("  miRNA names in the file resembling the families':", if (length(hits)) paste(head(hits, 10), collapse = ", ") else "none", "\n")
}

# ---- TarBase ---------------------------------------------------------------------------------------
cat(sprintf("\n== TarBase: %s ==\n", HUMAN_DB$tarbase))
if (!file.exists(HUMAN_DB$tarbase)) {
  cat("  FILE NOT FOUND -> step 1 skipped TarBase. Put the file here or set HUMAN_DB$tarbase in benchmark_config.R.\n")
} else {
  h <- file_header(HUMAN_DB$tarbase)
  cat("  columns:", paste(h, collapse = ", "), "\n")
  mcol <- pick_col(h, c("mirna_name", "mirna"), "TarBase miRNA"); gcol <- pick_col(h, "gene_id", "TarBase gene_id")
  tb <- data.table::fread(HUMAN_DB$tarbase, data.table = FALSE, select = c(mcol, gcol)); names(tb) <- c("mirna", "gene_id")
  stage("rows in file", nrow(tb))
  cat("  example miRNA names:", paste(head(unique(tb$mirna), 5), collapse = ", "), "| example gene_id:", paste(head(tb$gene_id, 3), collapse = ", "), "\n")
  tb <- filter(tb, mirna %in% members); stage("rows for the family members", nrow(tb))
  if (!nrow(tb)) similar_names(data.table::fread(HUMAN_DB$tarbase, data.table = FALSE, select = mcol)[[1]])
  tb$gene <- sub("\\..*$", "", tb$gene_id)
  stage("... with an Ensembl gene ID (ENSG)", sum(grepl("^ENSG", tb$gene)))
  stage("... whose gene is among the scored genes", sum(tb$gene %in% universe))
  stage("distinct scored target genes", n_distinct(tb$gene[tb$gene %in% universe]))
}

# ---- miRTarBase -------------------------------------------------------------------------------------
cat(sprintf("\n== miRTarBase: %s ==\n", HUMAN_DB$mirtarbase))
if (!file.exists(HUMAN_DB$mirtarbase)) {
  cat("  FILE NOT FOUND -> step 1 skipped miRTarBase. Put the file here or set HUMAN_DB$mirtarbase in benchmark_config.R.\n")
} else {
  mt <- readr::read_csv(HUMAN_DB$mirtarbase, show_col_types = FALSE)
  h <- names(mt)
  cat("  columns:", paste(h, collapse = ", "), "\n")
  if (ncol(mt) < 5) cat("  *** fewer than 5 columns: the file may not be comma-separated (e.g. tab- or semicolon-separated)\n")
  cols <- c(mirna = pick_col(h, "miRNA", "miRNA"), gene = pick_col(h, "Target Gene", "target gene"),
            species = pick_col(h, c("Species (Target Gene)", "Species..Target.Gene."), "target species"),
            support = pick_col(h, c("Support Type", "Support.Type"), "support type"))
  mt <- tibble(mirna = mt[[cols["mirna"]]], gene = mt[[cols["gene"]]], species = mt[[cols["species"]]], support = mt[[cols["support"]]])
  stage("rows in file", nrow(mt))
  cat("  species values (top):", paste(names(head(sort(table(mt$species), decreasing = TRUE), 4)), collapse = " | "), "\n")
  cat("  support types:", paste(names(table(mt$support)), collapse = " | "), "\n")
  cat("  example miRNA names:", paste(head(unique(mt$mirna), 5), collapse = ", "), "\n")
  x <- filter(mt, mirna %in% members); stage("rows for the family members", nrow(x))
  if (!nrow(x)) similar_names(mt$mirna)
  x <- filter(x, species == "Homo sapiens"); stage("... with target species exactly 'Homo sapiens'", nrow(x))
  x <- filter(x, support %in% c("Functional MTI", "Functional MTI (Weak)")); stage("... with a functional support type", nrow(x))
  x$gid <- symbol_to_gene(x$gene, utr, universe)
  stage("... whose gene symbol maps to a scored gene", sum(!is.na(x$gid)))
  stage("distinct scored target genes", n_distinct(na.omit(x$gid)))
}
cat("\nThe first stage that drops to 0 (or far below the one above it) is the cause.\n")
