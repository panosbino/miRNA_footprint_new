library(tidyverse)

# ---------------------------------------------------------------------------
# Exports mESC data for miTEA-HiRes, using the CUSTOM target-list bypass
# (same mechanism already built/verified for TarBase/TargetScan on human
# data) rather than the native miRTarBase pathway -- deliberately sidesteps
# needing to know the exact mouse miRTarBase filename or species-string
# convention, since the bypass works with any gene list regardless of
# species. Only requirement: mouse gene SYMBOL mapping, reused from
# mm.utr3.seqs.rds (same gid/gsym structure as the human version).
# ---------------------------------------------------------------------------

BASE_DIR <- "/cfs/klemming/projects/supr/naiss2024-6-235/miRNA_footprint_new"
TOOLS_DIR <- file.path(BASE_DIR, "tools")
MIREACT_DIR <- file.path(TOOLS_DIR, "miReact")
OUT_DIR <- file.path(BASE_DIR, "analysis/mESCs")
MITEA_INPUT_DIR <- file.path(OUT_DIR, "mitea_input")
dir.create(file.path(MITEA_INPUT_DIR, "counts_data"), recursive = TRUE, showWarnings = FALSE)

setwd(BASE_DIR)
pheno <- read.delim("./datasets/mESCs/phenotype.csv", sep = ",")
counts <- readRDS("datasets/mESCs/final_counts.rds")
counts <- as.matrix(counts[, colnames(counts) %in% pheno$cells])
cat(sprintf("Counts: %d genes x %d cells\n", nrow(counts), ncol(counts)))

# --- Mouse gene-symbol relabeling, same pattern as the HEK293 export -------
raw_seqs <- readRDS(file.path(MIREACT_DIR, "seqs", "mm.utr3.seqs.rds"))
gid_to_gsym <- raw_seqs %>%
  select(gid, gsym) %>%
  distinct() %>%
  filter(!is.na(gsym), gsym != "")
dup_gid <- sum(duplicated(gid_to_gsym$gid))
cat(sprintf("gid->gsym mapping: %d rows, %d gids with >1 candidate symbol (keeping first)\n",
            nrow(gid_to_gsym), dup_gid))
gid_to_gsym <- gid_to_gsym[!duplicated(gid_to_gsym$gid), ]

matched_idx <- match(rownames(counts), gid_to_gsym$gid)
n_mapped <- sum(!is.na(matched_idx))
cat(sprintf("Of %d genes, %d (%.1f%%) mapped to a gene symbol\n", nrow(counts), n_mapped, 100 * n_mapped / nrow(counts)))

counts_symbol <- counts[!is.na(matched_idx), ]
rownames(counts_symbol) <- gid_to_gsym$gsym[matched_idx[!is.na(matched_idx)]]
n_dup_symbols <- sum(duplicated(rownames(counts_symbol)))
cat(sprintf("%d duplicate gene-symbol rows (summing counts across duplicates)\n", n_dup_symbols))
counts_symbol_df <- as.data.frame(counts_symbol) %>%
  rownames_to_column("gene_symbol") %>%
  group_by(gene_symbol) %>%
  summarise(across(everything(), sum)) %>%
  column_to_rownames("gene_symbol")

out_file <- file.path(MITEA_INPUT_DIR, "counts_data", "counts_symbol.txt")
write.table(counts_symbol_df, out_file, sep = "\t", quote = FALSE, col.names = NA)
cat(sprintf("Wrote %s (%d genes x %d cells)\n", out_file, nrow(counts_symbol_df), ncol(counts_symbol_df)))

# --- Combined (3-family) and negative-control target lists, to symbols ----
targets_combined <- readRDS("resources/mouse/Targets_combined_top3_families_mESCs.rds")
targets_negctrl <- readRDS("resources/mouse/Targets_Negative_Control_mESCs.rds")

combined_symbols <- unique(na.omit(gid_to_gsym$gsym[match(targets_combined$ensembl_gene_id, gid_to_gsym$gid)]))
negctrl_symbols <- unique(na.omit(gid_to_gsym$gsym[match(targets_negctrl$ensembl_gene_id, gid_to_gsym$gid)]))

cat(sprintf("Combined target list: %d Ensembl IDs -> %d symbols\n", nrow(targets_combined), length(combined_symbols)))
cat(sprintf("Negative control list: %d Ensembl IDs -> %d symbols\n", nrow(targets_negctrl), length(negctrl_symbols)))

writeLines(combined_symbols, file.path(MITEA_INPUT_DIR, "combined_targets.txt"))
writeLines(negctrl_symbols, file.path(MITEA_INPUT_DIR, "negctrl_targets.txt"))
cat("Wrote combined_targets.txt and negctrl_targets.txt\n")
