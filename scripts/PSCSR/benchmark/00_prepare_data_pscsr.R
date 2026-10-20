# ==============================================================================
# PSCSR benchmark, step 0: prepare matrices and select the miRNAs to evaluate.
#
# mRNA (same recipe as the HEK scran matrix, build_annotated_scran_dataset.R):
#   gene IDs -> Ensembl (org.Hs.eg.db for symbols; duplicates summed),
#   protein-coding genes, genes with zeros in < MAX_ZERO_PROP of cells,
#   scran quickCluster + computeSumFactors, linear normalised counts.
# miRNA: counts normalised (MIRNA_NORM: median of ratios per cell line, or CPM);
#   names matched to TargetScan's
#   miR_Family_Info (human, 9606) to get the seed (nt 2-8) and family.
# Units: UNIT_LEVEL "family" -> one unit per seed; all human miRNAs with that
#   seed share its targets, and its measured expression is the sum of their
#   members' normalised counts. UNIT_LEVEL "mirna" -> one unit per miRNA.
# Selection: the top TOP_N_PER_TYPE UNITS by mean expression IN EACH CELL LINE,
#   separately; each line is later evaluated only on its own list (selection.tsv).
#   Ranking families rather than miRNAs keeps 10 distinct target sets per line
#   (K562's top-10 miRNAs are only 6 seed families).
#   Target-site motif for bayesReact / miReact = reverse complement of nt 2-8.
#   Cross-checked against the dataset's own annotation, whose "target_site"
#   column holds the SEED for 5p arms and the SITE for 3p arms (verified on
#   known miRNAs) -- so it is only used after that correction.
#
# Outputs (analysis/PSCSR/benchmark/data/):
#   scran_normalized_linear.rds  genes x cells (Ensembl)
#   unit_expression.rds          units x cells, measured miRNA expression (ground truth)
#   cells.tsv                    cell, cell_type
#   selection.tsv                per cell line: rank, unit, mean expression, measured members
#   units.tsv                    all selected units: seed, motif, member miRNAs (targets built once per unit)
#   mirnas_for_mirdip.txt        member miRNA names to query in mirDIP
#   midgd_counts.tsv             raw integer counts, all genes (Ensembl) -- miDGD input
#   mitea_input/counts_data/counts_symbol.txt                            -- miTEA input
# ==============================================================================

suppressPackageStartupMessages({ library(tidyverse); library(Matrix); library(scran); library(scuttle); library(org.Hs.eg.db) })
source(file.path(Sys.getenv("MIRNA_BASE_DIR", "~/Desktop/Projects/miRNA_footprint_new"),
                 "scripts/PSCSR/benchmark/config_pscsr.R"))
dir.create(file.path(DATA_DIR, "mitea_input/counts_data"), recursive = TRUE, showWarnings = FALSE)

# ---- Load and align cells ------------------------------------------------------------------
mrna  <- load_mrna_counts()
mirna <- load_mirna_counts()
meta  <- load_cell_meta()
cells <- Reduce(intersect, list(colnames(mrna), colnames(mirna), meta$cell))
cat(sprintf("Cells: mRNA %d, miRNA %d, metadata %d -> %d in all three\n",
            ncol(mrna), ncol(mirna), nrow(meta), length(cells)))
if (length(cells) < 0.5 * min(ncol(mrna), ncol(mirna))) stop("Fewer than half the cells match -- check cell naming across files.")
meta  <- meta[match(cells, meta$cell), ]
mrna  <- mrna[, cells]; mirna <- mirna[, cells]
print(table(meta$cell_type))

nonint <- mean(mrna != round(mrna))
cat(sprintf("mRNA: %d features; %.2f%% of entries non-integer%s\n", nrow(mrna), 100 * nonint,
            if (nonint > 0) " (multi-mapping weights; rounded for miDGD only)" else ""))

# ---- mRNA gene IDs -> Ensembl ----------------------------------------------------------------
ids <- sub("\\.\\d+$", "", rownames(mrna))
is_ens <- grepl("^ENSG\\d+$", ids)
cat(sprintf("Feature IDs: %d Ensembl, %d other (treated as gene symbols)\n", sum(is_ens), sum(!is_ens)))
if (any(!is_ens)) {
  valid <- intersect(unique(ids[!is_ens]), AnnotationDbi::keys(org.Hs.eg.db, "SYMBOL"))
  sym2ens <- if (length(valid)) suppressMessages(AnnotationDbi::mapIds(org.Hs.eg.db, keys = valid, keytype = "SYMBOL",
                                                                      column = "ENSEMBL", multiVals = "first")) else character(0)
  ids[!is_ens] <- unname(sym2ens[ids[!is_ens]])                    # non-symbols (spike-ins, transgenes) become NA and are dropped
  cat(sprintf("  symbols mapped to Ensembl: %d of %d\n", sum(!is.na(ids[!is_ens])), sum(!is_ens)))
}
keep <- !is.na(ids)
mrna <- rowsum(mrna[keep, , drop = FALSE], ids[keep])              # sums genes that map to the same Ensembl ID
cat(sprintf("mRNA after ID mapping: %d Ensembl genes\n", nrow(mrna)))

# miDGD input: all genes, raw integer counts
data.table::fwrite(data.table::as.data.table(round(mrna), keep.rownames = "gene_id"),
                   file.path(DATA_DIR, "midgd_counts.tsv"), sep = "\t")

# ---- scran matrix (HEK recipe) -----------------------------------------------------------------
pc <- readRDS(PROTEIN_CODING_FILE)$ensembl_gene_id
m_pc <- mrna[rownames(mrna) %in% pc, ]
m_pc <- m_pc[rowMeans(m_pc == 0) < MAX_ZERO_PROP, ]
cat(sprintf("Protein-coding genes with zeros in < %.0f%% of cells: %d\n", 100 * MAX_ZERO_PROP, nrow(m_pc)))
sce <- SingleCellExperiment::SingleCellExperiment(assays = list(counts = m_pc))
set.seed(SCRAN_SEED)
clusters <- quickCluster(sce)
cat("quickCluster vs cell type:\n"); print(table(cluster = clusters, type = meta$cell_type))
sce <- computeSumFactors(sce, clusters = clusters)
sce <- logNormCounts(sce, log = FALSE)
scran_norm <- as.matrix(SummarizedExperiment::assay(sce, "normcounts"))
saveRDS(scran_norm, file.path(DATA_DIR, "scran_normalized_linear.rds"))

# miTEA input: same gid -> symbol rule as export_for_mitea_rebuild.R
utr <- readRDS(file.path(MIREACT_DIR, "seqs/hs.utr3.seqs.rds"))
gid_to_gsym <- utr |> distinct(gid, gsym) |> filter(!is.na(gsym), gsym != "") |> filter(!duplicated(gid))
idx <- match(rownames(scran_norm), gid_to_gsym$gid)
cs <- rowsum(scran_norm[!is.na(idx), ], gid_to_gsym$gsym[idx[!is.na(idx)]])
write.table(as.data.frame(cs), file.path(DATA_DIR, "mitea_input/counts_data/counts_symbol.txt"),
            sep = "\t", quote = FALSE, col.names = NA)

# ---- miRNA: CPM and seed annotation ----------------------------------------------------------
rownames(mirna) <- standardise_mirna_names(rownames(mirna))
fam <- read_tsv(TS_FAMILY_FILE, show_col_types = FALSE)
fam <- fam |> filter(.data[[pick_col(fam, c("Species ID", "Species.ID"), "species")]] == 9606) |>
  transmute(mirbase_id = .data[[pick_col(fam, c("MiRBase ID", "miRBase ID"), "miRBase ID")]],
            ts_family  = .data[[pick_col(fam, c("miR family", "miR Family"), "family")]],
            seed       = toupper(gsub("U", "T", .data[[pick_col(fam, c("Seed+m8", "Seed + m8"), "seed")]])))
stopifnot(!anyDuplicated(fam$mirbase_id))

lib <- colSums(mirna)
if (MIRNA_NORM == "median_ratio") {
  stable <- rownames(mirna)[rowMeans(mirna > 0) >= 0.9]
  sf <- setNames(numeric(ncol(mirna)), colnames(mirna))
  for (ct in unique(meta$cell_type)) {                 # reference per cell line: line-level composition is not normalised away
    i <- meta$cell_type == ct
    L <- log(mirna[stable, i, drop = FALSE]); L[!is.finite(L)] <- NA
    ref <- rowMeans(L, na.rm = TRUE)
    sf[i] <- apply(L, 2, function(l) exp(median(l - ref, na.rm = TRUE)))
  }
  sf <- sf / exp(mean(log(sf)))
  cpm <- sweep(mirna, 2, sf * median(lib), "/") * 1e6   # scaled so values stay on a CPM-like scale
  cat(sprintf("miRNA normalisation: median of ratios over %d miRNAs detected in >= 90%% of cells; rho(size factor, library size) = %.2f\n",
              length(stable), cor(sf, lib, method = "spearman")))
} else {
  cpm <- sweep(mirna, 2, lib, "/") * 1e6
  cat("miRNA normalisation: CPM of the miRNA library\n")
}
matched <- rownames(cpm) %in% fam$mirbase_id
cat(sprintf("miRNA table: %d features, %d match a human miRBase ID in miR_Family_Info (%.1f%% of counts)\n",
            nrow(cpm), sum(matched), 100 * sum(mirna[matched, ]) / sum(mirna)))

# ---- Candidate units: every seed family (or miRNA) present in the table -------------------------------
fam_tab <- fam |> filter(mirbase_id %in% rownames(cpm)[matched])     # miRNAs measured here, with seed/family
if (UNIT_LEVEL == "family") {
  groups <- split(fam_tab$mirbase_id, fam_tab$seed)                 # seed -> measured members
  cand <- fam |> group_by(seed) |>
    summarise(ts_family = paste(unique(ts_family), collapse = "+"),
              members = paste(sort(unique(mirbase_id)), collapse = ";"), .groups = "drop") |>   # ALL human members (for databases)
    filter(seed %in% names(groups)) |>
    mutate(key = seed, unit = gsub("[/ ]", "_", ts_family))
} else {
  groups <- setNames(as.list(fam_tab$mirbase_id), fam_tab$mirbase_id)
  cand <- fam_tab |> transmute(seed, ts_family, members = mirbase_id, key = mirbase_id, unit = sub("^hsa-", "", mirbase_id))
}
stopifnot(!anyDuplicated(cand$unit))
expr_all <- t(vapply(groups[cand$key], function(m) colSums(cpm[m, , drop = FALSE]), numeric(ncol(cpm))))
dimnames(expr_all) <- list(cand$unit, colnames(cpm))

# ---- Top N units PER CELL LINE (each line is evaluated on its own list) ---------------------------------
sel <- map_dfr(CELL_TYPES, function(ct) {
  i <- meta$cell_type == ct
  all_top <- head(sort(rowMeans(cpm[, i, drop = FALSE]), decreasing = TRUE), TOP_N_PER_TYPE)
  miss <- setdiff(names(all_top), fam$mirbase_id)                    # naming problems stay visible
  if (length(miss)) cat(sprintf("  %s: top features without a miRBase/TargetScan match (skipped): %s\n", ct, paste(miss, collapse = ", ")))
  m <- head(sort(rowMeans(expr_all[, i, drop = FALSE]), decreasing = TRUE), TOP_N_PER_TYPE)
  tibble(cell_type = ct, rank = seq_along(m), unit = names(m), mean_expr = unname(m),
         frac_detected = unname(rowMeans(expr_all[names(m), i, drop = FALSE] > 0)),
         members_detected = vapply(groups[cand$key[match(names(m), cand$unit)]], paste, "", collapse = ";"))
})
units <- cand |> filter(unit %in% sel$unit) |> dplyr::select(unit, seed, ts_family, members)
units$motif <- vapply(units$seed, revcomp_dna, "")
units$members_detected <- vapply(groups[cand$key[match(units$unit, cand$unit)]], paste, "", collapse = ";")
unit_members <- setNames(strsplit(units$members, ";"), units$unit)
unit_expr <- expr_all[units$unit, , drop = FALSE]

# Cross-check motifs against the dataset annotation (target_site: seed for 5p, site for 3p)
if (file.exists(MIRNA_ANNO_FILE)) {
  ann <- read_tsv(MIRNA_ANNO_FILE, show_col_types = FALSE) |>
    mutate(site = ifelse(arm_origin == "5p", vapply(target_site, revcomp_dna, ""), toupper(target_site)))
  chk <- units |> mutate(mirna = strsplit(members_detected, ";")) |> unnest(mirna) |>
    mutate(ann_site = ann$site[match(mirna, ann$mature_mir_name)])
  bad <- chk |> filter(!is.na(ann_site), motif != ann_site)
  cat(sprintf("Motif cross-check vs %s: %d of %d measured members of the selected units agree%s\n", basename(MIRNA_ANNO_FILE),
              sum(chk$motif == chk$ann_site, na.rm = TRUE), sum(!is.na(chk$ann_site)), if (nrow(bad)) " -- MISMATCHES:" else ""))
  if (nrow(bad)) print(bad |> dplyr::select(unit, mirna, motif, ann_site))
}

cat(sprintf("\nTop %d %s units per cell line (%d distinct units overall):\n", TOP_N_PER_TYPE, UNIT_LEVEL, nrow(units)))
for (ct in CELL_TYPES) {
  s <- sel |> filter(cell_type == ct)
  cv <- apply(unit_expr[s$unit, meta$cell_type == ct, drop = FALSE], 1, \(v) sd(v) / mean(v))
  cat(sprintf("\n%s:\n", ct))
  print(s |> mutate(within_cv = round(unname(cv), 2), mean_expr = round(mean_expr)) |>
          dplyr::select(rank, unit, mean_expr, frac_detected, within_cv, members_detected), n = Inf, width = Inf)
}
cat("within_cv = coefficient of variation inside the line (low = little variation for activity to track)\n")

write_tsv(meta, file.path(DATA_DIR, "cells.tsv"))
write_tsv(sel, file.path(DATA_DIR, "selection.tsv"))
write_tsv(units, file.path(DATA_DIR, "units.tsv"))
saveRDS(unit_expr, file.path(DATA_DIR, "unit_expression.rds"))
writeLines(sort(unique(unlist(unit_members))), file.path(DATA_DIR, "mirnas_for_mirdip.txt"))
cat(sprintf("\nWrote prepared data to %s\nNext: query mirDIP for the miRNAs in mirnas_for_mirdip.txt and save the export as MIRDIP_FILE.\n", DATA_DIR))
