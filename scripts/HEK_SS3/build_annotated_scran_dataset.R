#!/usr/bin/env Rscript

# ==============================================================================
# Build annotated, scran-normalized single-cell dataset (NO Seurat)
#
# Merges:
#   - combine_plates_umap_UMIs.R  (load + combine zUMIs UMI counts across plates)
#   - annotate_umap_dox_gfp.R     (Dox concentration / induction group / GFP)
#   - compute_scran_matrix.R      (scran size-factor normalization)
#
# Seurat is removed entirely. Normalization is scran/scuttle-based throughout;
# dimensionality reduction and plotting use scater + ggplot2 directly on a
# SingleCellExperiment.
#
# MERGE NOTES / FIXES (on top of the ones already applied to the Seurat-based
# versions of the first two scripts):
#   1. Sparse matrices are kept as dgCMatrix throughout. The original
#      compute_scran_matrix.R did
#        readRDS(...) |> as.matrix() |> as.data.frame()
#      which densifies the whole count matrix -- expensive and unnecessary,
#      since SingleCellExperiment and scran both work natively with sparse
#      (dgCMatrix) assays.
#   2. The GFP gene is now resolved dynamically (reusing the multi-name /
#      ambiguity-reporting logic from the annotation script) instead of the
#      hardcoded `rownames(counts_raw) != "eGFP"` in the original scran
#      script. That hardcoded check would have silently done nothing if the
#      reporter gene were named e.g. "GFP" instead of "eGFP" -- it happened
#      to be harmless in practice because the later protein-coding filter
#      would have excluded it anyway, but that was a coincidence, not a
#      guarantee.
#   3. GFP is preserved from the RAW (unfiltered) counts and normalized
#      AFTER the fact using the same per-cell scran size factors used for
#      every other gene, so GFP is comparably normalized rather than absent
#      from the normalized dataset entirely (see note above the "GFP
#      normalization" section below).
#   4. Protein-coding gene list from biomaRt is cached locally after the
#      first query (previous script re-queried the public Ensembl server on
#      every run -- slow and fragile if Ensembl is briefly unreachable).
#   5. All previously-applied fixes are carried over: barcode/well
#      match-rate sanity checks, "Undefined"->"Unknown" recoding before
#      factor(), unified "No Dox"/"Dox Induced" label casing, zero-match
#      warnings, and a pre-plotting Dox/plate summary table.
#   6. No setwd() anywhere. The three original scripts each used a different,
#      machine-specific working directory (D:/microIMP,
#      /Volumes/Extreme SSD/microIMP, ~/Desktop/Projects/.../Manuscript_v0).
#      This version uses two explicit base-path variables (raw zUMIs data vs.
#      processed/output data) and builds every path with file.path(), so it
#      is portable across machines by editing two lines.
#   7. The Klemming HPC transfer reminder from compute_scran_matrix.R is
#      preserved verbatim at the end -- that's an operational note from your
#      colleague, not a bug, and easy to lose in a refactor.
# ==============================================================================

## ---- Config ---------------------------------------------------------------
# Where the raw per-plate zUMIs output lives (dgecounts.rds, barcode files,
# well maps) -- corresponds to "D:/microIMP" / plate subfolders in the
# original combine/annotate scripts.
RAW_DATA_DIR <- "/cfs/klemming/projects/snic/naiss2024-6-235/miRNA_footprint/datasets/HEK_SS3/raw"

# Where combined/normalized matrices and metadata are read from and written
# to -- corresponds to "../HEK_miR124_GFP/Data/" in compute_scran_matrix.R.
PROCESSED_DATA_DIR <- "/cfs/klemming/projects/snic/naiss2024-6-235/miRNA_footprint/datasets/HEK_SS3/processed"

# Where plots get written.
PLOT_DIR <- "/cfs/klemming/projects/snic/naiss2024-6-235/miRNA_footprint/analysis/HEK_SS3/qc/plots"

PLATES <- c("101", "103", "105")
SEQUENCING_PROJECT_PREFIX <- "P32156"

MIN_BARCODE_MATCH_FRAC <- 0.5
POSSIBLE_GFP_NAMES <- c("eGFP", "GFP", "Gfp", "EGFP", "EGfp")
MAX_ZERO_PROP <- 0.75          # from compute_scran_matrix.R
N_HVGS <- 2000                  # matches nfeatures used in the old Seurat step
N_PCA_DIMS <- 30
SEED <- 42
SCRAN_CLUSTER_SEED <- 123       # matches compute_scran_matrix.R's own seed

## ---- Libraries --------------------------------------------------------------
library(tidyverse)             # dplyr, stringr, ggplot2, etc. (as in the original scran script)
library(Matrix)
library(SingleCellExperiment)
library(scran)
library(scuttle)
library(scater)                 # PCA/UMAP + plotting helpers for SCE (replaces Seurat's RunPCA/RunUMAP)
library(biomaRt)
library(patchwork)
library(ggridges)

set.seed(SEED)

## ============================================================================
## 1. Load and combine raw UMI counts across plates (zUMIs output)
## ============================================================================

filter_count_data <- function(plate_id) {
  rds_file_path <- file.path(RAW_DATA_DIR,
                              paste0(SEQUENCING_PROJECT_PREFIX, "_", plate_id, ".dgecounts.rds"))
  barcode_file_path <- file.path(RAW_DATA_DIR,
                                  paste0("cell_barcodes_filtered_", plate_id, ".txt"))

  filtered_barcodes <- readLines(barcode_file_path)

  message(paste("Loading data for plate", plate_id))
  all_counts <- readRDS(rds_file_path)
  umi_counts_inex <- all_counts$umicount$inex$all   # dgCMatrix, kept sparse throughout

  retained_barcodes <- intersect(colnames(umi_counts_inex), filtered_barcodes)
  if (length(retained_barcodes) == 0) {
    stop(paste("No matching barcodes found for plate", plate_id))
  }
  match_frac <- length(retained_barcodes) / length(filtered_barcodes)
  if (match_frac < MIN_BARCODE_MATCH_FRAC) {
    warning(sprintf(
      "Plate %s: only %d / %d filtered barcodes (%.1f%%) matched the count matrix. Check barcode formatting.",
      plate_id, length(retained_barcodes), length(filtered_barcodes), 100 * match_frac
    ))
  }
  message(paste("Filtered", length(retained_barcodes), "barcodes for plate", plate_id))

  filtered_matrix <- umi_counts_inex[, retained_barcodes, drop = FALSE]
  colnames(filtered_matrix) <- paste0(plate_id, "_", colnames(filtered_matrix))
  return(filtered_matrix)
}

count_matrices <- lapply(PLATES, filter_count_data)
names(count_matrices) <- PLATES

all_genes <- lapply(count_matrices, rownames)
common_genes <- Reduce(intersect, all_genes)
message(paste("Found", length(common_genes), "common genes across all plates"))

count_matrices <- lapply(count_matrices, function(m) m[common_genes, , drop = FALSE])
combined_counts <- do.call(cbind, count_matrices)   # dgCMatrix, genes x cells

cell_ids <- colnames(combined_counts)
cell_metadata <- data.frame(
  cell_id = cell_ids,
  plate = sapply(strsplit(cell_ids, "_"), function(x) x[1]),
  row.names = cell_ids,
  stringsAsFactors = FALSE
)

# Preserve the same on-disk contract compute_scran_matrix.R originally relied
# on, in case any other script/collaborator still reads this file directly.
saveRDS(combined_counts, file.path(PROCESSED_DATA_DIR, "combined_counts_filtered_UMIs.rds"))
message("Saved combined_counts_filtered_UMIs.rds")

## ============================================================================
## 2. Dox concentration / induction group metadata
## ============================================================================

extract_dox_info <- function(plate_id) {
  barcode_file_path <- file.path(RAW_DATA_DIR,
                                  paste0("cell_barcodes_filtered_", plate_id, ".txt"))
  well_map_file_path <- file.path(RAW_DATA_DIR,
                                   paste0(SEQUENCING_PROJECT_PREFIX, "_", plate_id, ".well_barcodes.txt"))

  filtered_barcodes <- readLines(barcode_file_path)
  barcode_well_map <- read.delim(well_map_file_path, header = TRUE, sep = "\t")
  rownames(barcode_well_map) <- barcode_well_map$bc_set

  common_barcodes <- intersect(filtered_barcodes, rownames(barcode_well_map))
  if (length(common_barcodes) == 0) {
    warning(sprintf("Plate %s: NO barcodes matched the well map. Dox info will be entirely 'Unknown'.", plate_id))
  } else {
    match_frac <- length(common_barcodes) / length(filtered_barcodes)
    if (match_frac < MIN_BARCODE_MATCH_FRAC) {
      warning(sprintf(
        "Plate %s: only %d / %d filtered barcodes (%.1f%%) matched the well map.",
        plate_id, length(common_barcodes), length(filtered_barcodes), 100 * match_frac
      ))
    }
  }

  dox_df <- data.frame(
    barcode = common_barcodes,
    WellID = barcode_well_map[common_barcodes, "WellID"],
    row.names = common_barcodes
  )
  dox_df$Plate_Column <- as.numeric(gsub("^[A-Z]+0?(\\d+)$", "\\1", dox_df$WellID))

  # NOTE: plate 105 has a different column layout than 101/103; if a future
  # plate is added, revisit this if(plate_id == "105") check.
  if (plate_id == "105") {
    dox_df$Dox_Concentration <- case_when(
      dox_df$WellID %in% c("P23", "P24") ~ "Control",
      dox_df$Plate_Column >= 1  & dox_df$Plate_Column <= 5  ~ "0 Dox",
      dox_df$Plate_Column >= 6  & dox_df$Plate_Column <= 11 ~ "0.01 Dox",
      dox_df$Plate_Column >= 12 & dox_df$Plate_Column <= 17 ~ "0.1 Dox",
      dox_df$Plate_Column >= 18 & dox_df$Plate_Column <= 24 ~ "1 Dox",
      TRUE ~ "Undefined"
    )
  } else {
    dox_df$Dox_Concentration <- case_when(
      dox_df$WellID %in% c("P23", "P24") ~ "Control",
      dox_df$Plate_Column >= 1  & dox_df$Plate_Column <= 5  ~ "0 Dox",
      dox_df$Plate_Column >= 6  & dox_df$Plate_Column <= 11 ~ "0.01 Dox",
      dox_df$Plate_Column >= 12 & dox_df$Plate_Column <= 17 ~ "0.1 Dox",
      dox_df$Plate_Column >= 18 & dox_df$Plate_Column <= 23 ~ "1 Dox",
      dox_df$Plate_Column == 24 ~ "Control",
      TRUE ~ "Undefined"
    )
  }

  dox_df$cell_id <- paste0(plate_id, "_", dox_df$barcode)
  dox_df[, c("cell_id", "Dox_Concentration")]
}

message("Extracting Dox concentration info for each plate...")
all_dox_info <- do.call(rbind, lapply(PLATES, extract_dox_info))
rownames(all_dox_info) <- all_dox_info$cell_id

cell_metadata$Dox_Concentration <- "Unknown"
matched <- intersect(all_dox_info$cell_id, cell_metadata$cell_id)
cell_metadata[matched, "Dox_Concentration"] <- all_dox_info[matched, "Dox_Concentration"]

n_undefined <- sum(cell_metadata$Dox_Concentration == "Undefined")
if (n_undefined > 0) {
  warning(sprintf("%d cells had an 'Undefined' Dox_Concentration -- recoding to 'Unknown'.", n_undefined))
  cell_metadata$Dox_Concentration[cell_metadata$Dox_Concentration == "Undefined"] <- "Unknown"
}
cell_metadata$Dox_Concentration <- factor(
  cell_metadata$Dox_Concentration,
  levels = c("Control", "0 Dox", "0.01 Dox", "0.1 Dox", "1 Dox", "Unknown")
)

cell_metadata$Dox_Group <- "Unknown"
cell_metadata$Dox_Group[cell_metadata$Dox_Concentration %in% c("Control", "0 Dox")] <- "No Dox"
cell_metadata$Dox_Group[cell_metadata$Dox_Concentration %in% c("0.01 Dox", "0.1 Dox", "1 Dox")] <- "Dox Induced"
cell_metadata$Dox_Group <- factor(cell_metadata$Dox_Group, levels = c("No Dox", "Dox Induced", "Unknown"))

message("Dox_Concentration counts by plate:")
print(table(cell_metadata$plate, cell_metadata$Dox_Concentration))
n_unknown <- sum(cell_metadata$Dox_Concentration == "Unknown")
if (n_unknown > 0) {
  warning(sprintf("%d / %d cells (%.1f%%) have Unknown Dox_Concentration.",
                   n_unknown, nrow(cell_metadata), 100 * n_unknown / nrow(cell_metadata)))
}

## ============================================================================
## 3. Resolve the GFP gene ID (from the RAW, unfiltered counts)
## ============================================================================

egfp_gene_id <- NULL
for (name in POSSIBLE_GFP_NAMES) {
  if (name %in% rownames(combined_counts)) {
    egfp_gene_id <- name
    message(paste("Found GFP gene with exact ID:", egfp_gene_id))
    break
  }
}
if (is.null(egfp_gene_id)) {
  gfp_matches <- grep("GFP|gfp|Gfp", rownames(combined_counts), value = TRUE)
  if (length(gfp_matches) == 0) stop("No GFP gene found in the dataset!")
  if (length(gfp_matches) > 1) {
    warning(sprintf("Multiple GFP-like gene names found: %s. Using '%s'.",
                     paste(gfp_matches, collapse = ", "), gfp_matches[1]))
  }
  egfp_gene_id <- gfp_matches[1]
}
gfp_raw_counts <- setNames(as.numeric(combined_counts[egfp_gene_id, cell_metadata$cell_id]), cell_metadata$cell_id)

## ============================================================================
## 4. Prepare matrix for scran normalization (protein-coding, non-spike genes)
## ============================================================================

counts_raw <- combined_counts
rownames(counts_raw) <- str_split_i(rownames(counts_raw), pattern = "\\.", i = 1)  # strip Ensembl version suffix

counts_raw <- counts_raw[rownames(counts_raw) != egfp_gene_id, , drop = FALSE]
counts_raw <- counts_raw[!grepl("^ERCC", rownames(counts_raw)), , drop = FALSE]

pc_cache_path <- file.path(PROCESSED_DATA_DIR, "biomart_protein_coding_ids.rds")
if (file.exists(pc_cache_path)) {
  message("Loading cached protein-coding gene ID list...")
  pc_query <- readRDS(pc_cache_path)
} else {
  message("Querying Ensembl biomaRt for protein-coding gene IDs (first run only)...")
  ensembl_hs <- useEnsembl(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")
  pc_query <- getBM(attributes = c("ensembl_gene_id"),
                     filters = "biotype", values = "protein_coding", mart = ensembl_hs)
  saveRDS(pc_query, pc_cache_path)
}

counts_raw_PC <- counts_raw[rownames(counts_raw) %in% pc_query$ensembl_gene_id, , drop = FALSE]
counts_raw_PC <- counts_raw_PC[rowSums(is.na(counts_raw_PC)) == 0, , drop = FALSE]
cat(sprintf("Protein-coding genes: %d\n", nrow(counts_raw_PC)))

zero_prop <- Matrix::rowMeans(counts_raw_PC == 0)
counts_raw_PC_filt <- counts_raw_PC[zero_prop < MAX_ZERO_PROP, , drop = FALSE]
cat(sprintf("After zero-proportion filter: %d genes\n", nrow(counts_raw_PC_filt)))

## ============================================================================
## 5. scran normalization
## ============================================================================

sce <- SingleCellExperiment(assays = list(counts = counts_raw_PC_filt),
                             colData = cell_metadata[colnames(counts_raw_PC_filt), ])

set.seed(SCRAN_CLUSTER_SEED)
clusters <- quickCluster(sce)
sce <- computeSumFactors(sce, clusters = clusters)
sce <- logNormCounts(sce, log = FALSE)   # linear scale, matching the original script's design choice
scran_norm <- as.matrix(assay(sce, "normcounts"))

saveRDS(scran_norm, file.path(PROCESSED_DATA_DIR, "scran_normalized_linear.rds"))
cat(sprintf("Saved scran_normalized_linear.rds: %d genes x %d cells\n", nrow(scran_norm), ncol(scran_norm)))
cat("\n*** REMINDER: this file needs to be transferred to Klemming's ./Data/ directory\n",
    "before any Klemming-based script (miReact/bayesReact/miTEA comparisons) can use it. ***\n")

## ---- GFP normalization -----------------------------------------------------
# GFP was excluded from the matrix used to COMPUTE size factors (correct --
# a transgene shouldn't influence normalization), but the resulting
# scran_normalized_linear.rds therefore contains no GFP values at all. Since
# scran's per-cell size factors can be applied to any gene, apply them here
# to the raw GFP counts so a comparably-normalized GFP value is still
# available for the plots below. This step didn't exist in either original
# script -- flagging it as new, not a restored behavior.
size_factors <- sizeFactors(sce)
gfp_normalized <- gfp_raw_counts[colnames(sce)] / size_factors
cell_metadata$GFP_normalized <- NA_real_
cell_metadata[colnames(sce), "GFP_normalized"] <- gfp_normalized

# Also persist the raw GFP counts (indexed by ALL cells, not just those
# retained in sce, i.e. same set as combined_counts) so downstream scripts
# (e.g. a GFP-fluorescence correlation QC script) can compare raw vs.
# normalized GFP against FACS data without needing to reload combined_counts.
cell_metadata$GFP_raw_counts <- NA_real_
cell_metadata[names(gfp_raw_counts), "GFP_raw_counts"] <- gfp_raw_counts

## ============================================================================
## 6. Dimensionality reduction (QC / batch-effect check only)
## ============================================================================
# As agreed: this embedding is a quick sanity check (batch effects, Dox,
# GFP), not the "official" quantification -- that stays on the linear
# scran_normalized_linear.rds saved above. log1p() here is local/temporary,
# used only to make HVG selection and PCA well-behaved (scran's variance
# modeling expects roughly log-scale data); it is not written back to the
# saved normalized matrix.
logcounts(sce) <- log1p(assay(sce, "normcounts"))

dec <- modelGeneVar(sce)
hvgs <- getTopHVGs(dec, n = N_HVGS)
sce_hvg <- sce[hvgs, ]

set.seed(SEED)
sce_hvg <- runPCA(sce_hvg, ncomponents = N_PCA_DIMS)
sce_hvg <- runUMAP(sce_hvg, dimred = "PCA")

reducedDim(sce, "PCA") <- reducedDim(sce_hvg, "PCA")
reducedDim(sce, "UMAP") <- reducedDim(sce_hvg, "UMAP")

umap_df <- as.data.frame(reducedDim(sce, "UMAP"))
colnames(umap_df) <- c("UMAP1", "UMAP2")
umap_df$cell_id <- colnames(sce)
plot_df <- left_join(umap_df, cell_metadata, by = "cell_id")

# Persist UMAP coordinates into cell_metadata itself (not just the local
# plot_df) so anything saved to disk carries the embedding along.
cell_metadata$UMAP1 <- NA_real_
cell_metadata$UMAP2 <- NA_real_
cell_metadata[umap_df$cell_id, c("UMAP1", "UMAP2")] <- umap_df[, c("UMAP1", "UMAP2")]

## ============================================================================
## 7. Plots
## ============================================================================

theme_set(theme_minimal())

# -- Batch check: colored by plate --
p_batch <- ggplot(plot_df, aes(UMAP1, UMAP2, color = plate)) +
  geom_point(size = 2, alpha = 0.8) +
  scale_color_brewer(palette = "Dark2") +
  ggtitle("UMAP visualization of UMI reads colored by plate")
ggsave(file.path(PLOT_DIR, "umap_plate_batches_UMIs.png"), p_batch, width = 10, height = 8, dpi = 300)

# -- Dox concentration --
p_dox <- ggplot(plot_df, aes(UMAP1, UMAP2, color = Dox_Concentration)) +
  geom_point(size = 2, alpha = 0.8) +
  scale_color_brewer(palette = "Dark2") +
  ggtitle("UMAP visualization colored by Dox concentration")
ggsave(file.path(PLOT_DIR, "umap_dox_concentration.png"), p_dox, width = 10, height = 8, dpi = 300)

p_dox_facet <- ggplot(plot_df, aes(UMAP1, UMAP2, color = Dox_Concentration)) +
  geom_point(size = 1.5, alpha = 0.8) +
  facet_wrap(~ Dox_Concentration) +
  scale_color_brewer(palette = "Dark2") +
  ggtitle("UMAP visualization faceted by Dox concentration")
ggsave(file.path(PLOT_DIR, "umap_dox_concentration_faceted.png"), p_dox_facet, width = 16, height = 8, dpi = 300)

p_plate_dox <- (p_batch + ggtitle("By Plate")) | (p_dox + ggtitle("By Dox Concentration"))
ggsave(file.path(PLOT_DIR, "umap_plate_and_dox.png"), p_plate_dox, width = 16, height = 8, dpi = 300)

# -- Binary Dox induction group --
p_dox_binary <- ggplot(plot_df, aes(UMAP1, UMAP2, color = Dox_Group)) +
  geom_point(size = 2, alpha = 0.8) +
  scale_color_manual(values = c("No Dox" = "blue", "Dox Induced" = "red", "Unknown" = "gray")) +
  ggtitle("UMAP visualization by Dox induction status")
ggsave(file.path(PLOT_DIR, "umap_dox_binary.png"), p_dox_binary, width = 10, height = 8, dpi = 300)

p_dox_binary_split <- ggplot(plot_df, aes(UMAP1, UMAP2, color = Dox_Group)) +
  geom_point(size = 1.5, alpha = 0.8) +
  facet_wrap(~ plate) +
  scale_color_manual(values = c("No Dox" = "blue", "Dox Induced" = "red", "Unknown" = "gray")) +
  ggtitle("UMAP visualization by Dox induction status, split by plate")
ggsave(file.path(PLOT_DIR, "umap_dox_binary_by_plate.png"), p_dox_binary_split, width = 16, height = 8, dpi = 300)

# -- GFP expression --
percentile_99 <- quantile(plot_df$GFP_normalized, 0.99, na.rm = TRUE)
percentile_95 <- quantile(plot_df$GFP_normalized, 0.95, na.rm = TRUE)
message("GFP (normalized) expression summary:")
message(paste("  Max value:", max(plot_df$GFP_normalized, na.rm = TRUE)))
message(paste("  99th percentile:", percentile_99))
message(paste("  95th percentile:", percentile_95))

plot_gfp_umap <- function(cutoff, title, filename) {
  p <- ggplot(plot_df, aes(UMAP1, UMAP2, color = pmin(GFP_normalized, cutoff))) +
    geom_point(size = 2, alpha = 0.8) +
    scale_color_distiller(palette = "YlOrRd", direction = 1, name = "GFP\n(normalized)") +
    ggtitle(title)
  ggsave(file.path(PLOT_DIR, filename), p, width = 10, height = 8, dpi = 300)
  p
}
gfp_umap_99 <- plot_gfp_umap(percentile_99, paste(egfp_gene_id, "expression (99th percentile cutoff)"),
                              "umap_gfp_expression_99percentile.png")
gfp_umap_95 <- plot_gfp_umap(percentile_95, paste(egfp_gene_id, "expression (95th percentile cutoff)"),
                              "umap_gfp_expression_95percentile.png")

p_gfp_by_dox <- ggplot(plot_df %>% filter(Dox_Concentration != "Unknown"),
                        aes(UMAP1, UMAP2, color = pmin(GFP_normalized, percentile_99))) +
  geom_point(size = 1.5, alpha = 0.8) +
  facet_wrap(~ Dox_Concentration, nrow = 1) +
  scale_color_distiller(palette = "YlOrRd", direction = 1, name = "GFP\n(normalized)") +
  ggtitle(paste("UMAP colored by", egfp_gene_id, "expression by Dox concentration"))
ggsave(file.path(PLOT_DIR, "umap_gfp_expression_by_dox_improved.png"), p_gfp_by_dox, width = 24, height = 8, dpi = 300)

p_gfp_by_binary <- ggplot(plot_df %>% filter(Dox_Group != "Unknown"),
                           aes(UMAP1, UMAP2, color = pmin(GFP_normalized, percentile_99))) +
  geom_point(size = 2, alpha = 0.8) +
  facet_wrap(~ Dox_Group) +
  scale_color_distiller(palette = "YlOrRd", direction = 1, name = "GFP\n(normalized)") +
  ggtitle(paste("UMAP colored by", egfp_gene_id, "expression by Dox induction"))
ggsave(file.path(PLOT_DIR, "umap_gfp_expression_by_binary_improved.png"), p_gfp_by_binary, width = 18, height = 8, dpi = 300)

gfp_violin <- ggplot(plot_df, aes(x = Dox_Concentration, y = pmin(GFP_normalized, percentile_99 * 1.1),
                                   fill = Dox_Concentration)) +
  geom_violin(scale = "width") +
  scale_fill_brewer(palette = "Dark2") +
  labs(y = paste(egfp_gene_id, "expression (normalized)"),
       title = paste(egfp_gene_id, "expression by Dox concentration"))
ggsave(file.path(PLOT_DIR, "violin_gfp_expression_by_dox.png"), gfp_violin, width = 10, height = 6, dpi = 300)

gfp_ridge <- ggplot(plot_df, aes(x = pmin(GFP_normalized, percentile_99 * 1.1),
                                  y = Dox_Concentration, fill = Dox_Concentration)) +
  geom_density_ridges(alpha = 0.7) +
  scale_fill_brewer(palette = "Dark2") +
  labs(x = paste(egfp_gene_id, "Expression (normalized)"), y = "Dox Concentration",
       title = paste(egfp_gene_id, "expression distribution by Dox concentration"))
ggsave(file.path(PLOT_DIR, "ridge_gfp_expression_by_dox.png"), gfp_ridge, width = 10, height = 6, dpi = 300)

## ============================================================================
## 8. Save final objects
## ============================================================================
# The SCE (with PCA/UMAP embeddings, scran-normalized counts, and Dox/GFP
# colData) is now the single object downstream scripts should load --
# replaces the Seurat object entirely.
saveRDS(sce, file.path(PROCESSED_DATA_DIR, "sce_annotated.rds"))
saveRDS(cell_metadata, file.path(PROCESSED_DATA_DIR, "cell_metadata.rds"))
write.csv(cell_metadata, file.path(PROCESSED_DATA_DIR, "cell_metadata.csv"), row.names = FALSE)

message("Analysis complete. scran-normalized matrix, annotated SCE, cell metadata, and QC plots saved to: ",
        PROCESSED_DATA_DIR)
