# ==============================================================================
# PSCSR-seq V2, step 0: prepare matrices and choose the miRNA families to test.
#
# mRNA  : Ensembl IDs, protein-coding genes detected in > 25% of cells, scran
#         size factors (quickCluster + computeSumFactors), linear scale --
#         the recipe of the HEK scran matrix.
# miRNA : normalised by median of ratios per cell line (or CPM, PSCSR_MIRNA_NORM).
#         CPM is compositional: the miR-17 family is ~1/3 of K562 reads.
# Families: miRNAs sharing a seed (nt 2-8, from TargetScan's miR_Family_Info)
#         share their targets, so one family is one test. Its measured
#         expression = sum of its detected members.
# Selection: the PSCSR_TOP_N most expressed families IN EACH cell line; each
#         line is evaluated only on its own list.
#
# Outputs (analysis/PSCSR/benchmark_v3/data/):
#   scran_normalized_linear.rds  genes x cells
#   family_expression.rds        families x cells (ground truth)
#   cells.tsv, families.tsv, selection.tsv
#   mirnas_for_mirdip.txt        all members of the selected families AND of the control family
#   midgd_counts.tsv             raw integer counts, all genes (miDGD input)
#   mitea_input/counts_data/counts_symbol.txt
# ==============================================================================
source(file.path(Sys.getenv("MIRNA_BASE_DIR", "/cfs/klemming/projects/supr/naiss2024-6-235/miRNA_footprint_new"),
                 "scripts/benchmark_v3/benchmark_config.R"))
source(file.path(BASE_DIR, "scripts/benchmark_v3/benchmark_lib.R"))
suppressPackageStartupMessages({ library(scran); library(scuttle); library(org.Hs.eg.db) })
ds <- get_dataset("PSCSR")
data_dir <- file.path(ds$out_dir, "data")
dir.create(file.path(data_dir, "mitea_input/counts_data"), recursive = TRUE, showWarnings = FALSE)
inp <- function(f) file.path(PSCSR_INPUT$dir, PSCSR_INPUT[[f]])

# ---- 1. Load and align cells ----------------------------------------------------------------
mrna  <- read_count_matrix(inp("mrna_counts"))
mirna <- read_count_matrix(inp("mirna_counts"))
anno  <- data.table::fread(inp("cell_anno"), data.table = FALSE)
cells <- tibble(cell = as.character(anno[[PSCSR_INPUT$cell_col]]), cell_line = as.character(anno[[PSCSR_INPUT$line_col]]))
cells <- filter(cells, cell %in% colnames(mrna), cell %in% colnames(mirna))
stopifnot(all(cells$cell_line %in% PSCSR_CELL_LINES))
mrna <- mrna[, cells$cell]; mirna <- mirna[, cells$cell]
cat(sprintf("Cells with mRNA, miRNA and a cell line: %d\n", nrow(cells))); print(table(cells$cell_line))

# ---- 2. mRNA: gene IDs, miDGD input ----------------------------------------------------------------
rownames(mrna) <- sub("\\.\\d+$", "", rownames(mrna))
if (!all(grepl("^ENSG", rownames(mrna)))) stop("Expected Ensembl gene IDs in the mRNA table")
cat(sprintf("mRNA: %d genes; %.2f%% of entries non-integer (multi-mapping weights; rounded only for miDGD)\n",
            nrow(mrna), 100 * mean(mrna != round(mrna))))
data.table::fwrite(data.table::as.data.table(round(mrna), keep.rownames = "gene_id"), file.path(data_dir, "midgd_counts.tsv"), sep = "\t")

# ---- 3. mRNA: scran normalisation (HEK recipe) -------------------------------------------------------
protein_coding <- readRDS(PSCSR_PROTEIN_CODING)$ensembl_gene_id
m <- mrna[rownames(mrna) %in% protein_coding, ]
m <- m[rowMeans(m == 0) < COEXP_MAX_ZERO_FRAC, ]
cat(sprintf("Protein-coding genes detected in > %.0f%% of cells: %d\n", 100 * (1 - COEXP_MAX_ZERO_FRAC), nrow(m)))
sce <- SingleCellExperiment::SingleCellExperiment(assays = list(counts = m))
set.seed(123)                                           # as build_annotated_scran_dataset.R
clusters <- quickCluster(sce)
sce <- logNormCounts(computeSumFactors(sce, clusters = clusters), log = FALSE)
expr <- as.matrix(SummarizedExperiment::assay(sce, "normcounts"))
saveRDS(expr, file.path(data_dir, "scran_normalized_linear.rds"))

utr <- load_utr_table("hs")                             # miTEA input: gene -> symbol as the other miTEA exports
g2s <- utr |> distinct(gid, gsym) |> filter(!is.na(gsym), gsym != "") |> filter(!duplicated(gid))
idx <- match(rownames(expr), g2s$gid)
write.table(as.data.frame(rowsum(expr[!is.na(idx), ], g2s$gsym[idx[!is.na(idx)]])),
            file.path(data_dir, "mitea_input/counts_data/counts_symbol.txt"), sep = "\t", quote = FALSE, col.names = NA)

# ---- 4. miRNA normalisation ----------------------------------------------------------------------------------
library_size <- colSums(mirna)
if (PSCSR_MIRNA_NORM == "median_ratio") {
  stable <- rownames(mirna)[rowMeans(mirna > 0) >= 0.9]          # miRNAs detected in nearly every cell
  size_factor <- setNames(numeric(ncol(mirna)), colnames(mirna))
  for (line in PSCSR_CELL_LINES) {                                 # reference per line: line-level differences are kept
    i <- cells$cell_line == line
    logc <- log(mirna[stable, i, drop = FALSE]); logc[!is.finite(logc)] <- NA
    reference <- rowMeans(logc, na.rm = TRUE)
    size_factor[i] <- apply(logc, 2, \(l) exp(median(l - reference, na.rm = TRUE)))
  }
  size_factor <- size_factor / exp(mean(log(size_factor)))
  mirna_norm <- sweep(mirna, 2, size_factor * median(library_size), "/") * 1e6
  cat(sprintf("miRNA: median of ratios over %d miRNAs detected in >= 90%% of cells\n", length(stable)))
} else {
  mirna_norm <- sweep(mirna, 2, library_size, "/") * 1e6
  cat("miRNA: CPM\n")
}

# ---- 5. miRNA families ---------------------------------------------------------------------------------------
family_info <- read_family_info()
measured <- intersect(rownames(mirna_norm), family_info$mirbase_id)
cat(sprintf("miRNA table: %d miRNAs, %d with a seed in miR_Family_Info (%.1f%% of counts)\n", nrow(mirna), length(measured),
            100 * sum(mirna[measured, ]) / sum(mirna)))
seed_of_measured <- family_info$seed[match(measured, family_info$mirbase_id)]
family_table <- family_info |> filter(seed %in% seed_of_measured) |> group_by(seed) |>
  summarise(family = gsub("[/ ]", "_", paste(unique(ts_family), collapse = "+")),
            members = list(sort(unique(mirbase_id))), .groups = "drop") |>
  mutate(members_detected = map(members, \(x) intersect(x, measured)), site = vapply(seed, revcomp_dna, ""))
family_expr <- t(vapply(family_table$members_detected, \(x) colSums(mirna_norm[x, , drop = FALSE]), numeric(ncol(mirna_norm))))
dimnames(family_expr) <- list(family_table$family, colnames(mirna_norm))

# ---- 6. Top families per cell line -------------------------------------------------------------------------------
selection <- map_dfr(PSCSR_CELL_LINES, function(line) {
  i <- cells$cell_line == line
  top <- head(sort(rowMeans(family_expr[, i, drop = FALSE]), decreasing = TRUE), PSCSR_TOP_N)
  tibble(cell_line = line, rank = seq_along(top), family = names(top), mean_expression = unname(top),
         within_line_cv = apply(family_expr[names(top), i, drop = FALSE], 1, \(v) sd(v) / mean(v)))
})
families <- family_table |> filter(family %in% selection$family)

# Cross-check sites against the dataset's annotation (its target_site column is
# the SEED for 5p arms and the SITE for 3p arms)
ann <- read_tsv(inp("mirna_anno"), show_col_types = FALSE) |>
  mutate(site = ifelse(arm_origin == "5p", vapply(target_site, revcomp_dna, ""), toupper(target_site)))
chk <- families |> dplyr::select(family, site, members_detected) |> unnest(members_detected) |>
  mutate(annotated = ann$site[match(members_detected, ann$mature_mir_name)]) |> filter(!is.na(annotated))
cat(sprintf("Site cross-check vs the dataset annotation: %d of %d agree\n", sum(chk$site == chk$annotated), nrow(chk)))
if (any(chk$site != chk$annotated)) print(filter(chk, site != annotated))

for (line in PSCSR_CELL_LINES) {
  cat(sprintf("\n%s:\n", line))
  print(selection |> filter(cell_line == line) |> dplyr::select(-cell_line) |> mutate(mean_expression = round(mean_expression), within_line_cv = round(within_line_cv, 2)),
        n = Inf)
}

# ---- 7. Write -------------------------------------------------------------------------------------------------------
control <- human_family_of(ds$control_mirna, family_info)
write_tsv(cells, file.path(data_dir, "cells.tsv"))
write_tsv(selection, file.path(data_dir, "selection.tsv"))
write_tsv(families |> mutate(members = map_chr(members, paste, collapse = ";"),
                             members_detected = map_chr(members_detected, paste, collapse = ";")),
          file.path(data_dir, "families.tsv"))
saveRDS(family_expr[families$family, , drop = FALSE], file.path(data_dir, "family_expression.rds"))
writeLines(sort(unique(c(unlist(families$members), control$members))), file.path(data_dir, "mirnas_for_mirdip.txt"))
cat(sprintf("\n%d families selected (%d per line). Query mirDIP for the %d miRNAs in mirnas_for_mirdip.txt (includes the control %s family).\n",
            nrow(families), PSCSR_TOP_N, length(unique(c(unlist(families$members), control$members))), ds$control_mirna))
