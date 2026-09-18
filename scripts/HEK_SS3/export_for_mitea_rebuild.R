library(tidyverse)

# ---------------------------------------------------------------------------
# Rebuilt under miRNA_footprint_new, using scran-normalized data. Conda env
# (mitea, Python 3.8) is confirmed intact -- no Python-side setup needed.
# All known fixes from the original build reapplied (see inline notes).
# ---------------------------------------------------------------------------

BASE_DIR <- "/cfs/klemming/projects/supr/naiss2024-6-235/miRNA_footprint_new"
TOOLS_DIR <- file.path(BASE_DIR, "tools")
MIREACT_DIR <- file.path(TOOLS_DIR, "miReact")
PROCESSED_DIR <- file.path(BASE_DIR, "datasets/HEK_SS3/processed")
OUT_DIR <- file.path(BASE_DIR, "analysis/HEK_SS3/comparisons")
MITEA_INPUT_DIR <- file.path(OUT_DIR, "mitea_input")
dir.create(file.path(MITEA_INPUT_DIR, "counts_data"), recursive = TRUE, showWarnings = FALSE)

source(file.path(BASE_DIR, "scripts/Utils.R")) 

# --- Scran-normalized counts -------------------------------------------------
counts <- readRDS(file.path(PROCESSED_DIR, "scran_normalized_linear.rds"))
counts <- as.matrix(counts)
cat(sprintf("Counts: %d genes x %d cells (scran-normalized)\n", nrow(counts), ncol(counts)))

# --- Gene-symbol relabeling, reusing miReact's bundled mapping (unchanged logic) ---
raw_seqs <- readRDS(file.path(MIREACT_DIR, "seqs", "hs.utr3.seqs.rds"))
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
cat(sprintf("Of %d genes in counts matrix, %d (%.1f%%) mapped to a gene symbol\n",
            nrow(counts), n_mapped, 100 * n_mapped / nrow(counts)))

counts_symbol <- counts[!is.na(matched_idx), ]
rownames(counts_symbol) <- gid_to_gsym$gsym[matched_idx[!is.na(matched_idx)]]

n_dup_symbols <- sum(duplicated(rownames(counts_symbol)))
cat(sprintf("%d duplicate gene-symbol rows after mapping (summing counts across duplicates)\n", n_dup_symbols))
counts_symbol_df <- as.data.frame(counts_symbol) %>%
  rownames_to_column("gene_symbol") %>%
  group_by(gene_symbol) %>%
  summarise(across(everything(), sum)) %>%
  column_to_rownames("gene_symbol")

# Isolated in its OWN subdirectory (fix reapplied): process_data() globs
# for *.txt files in its data_path directory and errors if it finds more
# than one -- must not share a folder with tarbase_124_targets.txt below.
out_file <- file.path(MITEA_INPUT_DIR, "counts_data", "counts_symbol.txt")
write.table(counts_symbol_df, out_file, sep = "\t", quote = FALSE, col.names = NA)
cat(sprintf("Wrote %s (%d genes x %d cells)\n", out_file, nrow(counts_symbol_df), ncol(counts_symbol_df)))

# --- TarBase targets, for the TarBase-matched comparison arm ---------------
tar <- readRDS(file.path(MIREACT_DIR, "data", "tarbase.rds"))
tar_124 <- tar[tar$mirna == "hsa-miR-124-3p" & tar$species == "Homo sapiens" &
                 tar$up_down == "DOWN" & !is.na(tar$up_down), ]
tarbase_symbols <- unique(sub("\\(hsa\\)$", "", tar_124$geneName))
cat(sprintf("TarBase hsa-miR-124-3p DOWN targets: %d unique gene symbols\n", length(tarbase_symbols)))

writeLines(tarbase_symbols, file.path(MITEA_INPUT_DIR, "tarbase_124_targets.txt"))
cat(sprintf("Wrote %s (%d symbols)\n", file.path(MITEA_INPUT_DIR, "tarbase_124_targets.txt"), length(tarbase_symbols)))

# --- TargetScan targets, both full candidate list and top-200 --------------
# Reuses gid_to_gsym (already built above for the counts relabeling) --
# same mapping table, not a separately-derived one, for consistency with
# every other identifier-mapping step in this project.
N_TOP_TARGETS <- 200
targetscan_all <- readRDS(file.path(BASE_DIR, "resources/Targets__combined_124-3p_124-3p.2_506-3p.rds"))
targetscan_sorted <- targetscan_all[order(targetscan_all$Cumulative.weighted.context...score), ]
targets_top200_gid <- targetscan_sorted[1:N_TOP_TARGETS, ]$ensembl_gene_id

targetscan_all_symbols <- unique(na.omit(gid_to_gsym$gsym[match(targetscan_all$ensembl_gene_id, gid_to_gsym$gid)]))
targetscan_top200_symbols <- unique(na.omit(gid_to_gsym$gsym[match(targets_top200_gid, gid_to_gsym$gid)]))

cat(sprintf("TargetScan ALL candidates: %d Ensembl IDs -> %d unique gene symbols (%d lost in mapping)\n",
            nrow(targetscan_all), length(targetscan_all_symbols), nrow(targetscan_all) - length(targetscan_all_symbols)))
cat(sprintf("TargetScan TOP %d: %d Ensembl IDs -> %d unique gene symbols (%d lost in mapping)\n",
            N_TOP_TARGETS, length(targets_top200_gid), length(targetscan_top200_symbols),
            length(targets_top200_gid) - length(targetscan_top200_symbols)))

writeLines(targetscan_all_symbols, file.path(MITEA_INPUT_DIR, "targetscan_all_124_targets.txt"))
writeLines(targetscan_top200_symbols, file.path(MITEA_INPUT_DIR, "targetscan_top200_124_targets.txt"))
cat("Wrote targetscan_all_124_targets.txt and targetscan_top200_124_targets.txt\n")
