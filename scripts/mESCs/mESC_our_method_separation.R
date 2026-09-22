library(tidyverse)

BASE_DIR <- "/cfs/klemming/projects/supr/naiss2024-6-235/miRNA_footprint_new"
OUT_DIR <- file.path(BASE_DIR, "analysis/mESCs")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

source(file.path(BASE_DIR, "scripts/mESCs/mESC_separation_utils.R"))

# --- Data (exactly matching separate_KO_control.R's loading) ---------------
setwd(BASE_DIR)
pheno <- read.delim("./datasets/mESCs/phenotype.csv", sep = ",")
pheno$Exp[pheno$Exp == "WT"] <- "Control"
kos <- pheno[pheno$Exp == "KO", ]$cells

counts <- readRDS("datasets/mESCs/final_counts.rds")
counts <- counts[, colnames(counts) %in% pheno$cells]

targets_all <- readRDS("resources/mouse/Targets_combined_top3_families_mESCs.rds")
targets_all <- targets_all %>%
  group_by(ensembl_gene_id) %>%
  mutate(score = sum(Cumulative.weighted.context...score)) %>%
  ungroup() %>%
  distinct(pick(ensembl_gene_id), .keep_all = TRUE)

cat(sprintf("Combined top-3-families target list: %d unique genes\n", nrow(targets_all)))

targets_top200 <- targets_all %>% arrange(score) %>% slice_head(n = 200)  # most negative/strongest score first, matching TargetScan convention used throughout this project
cat(sprintf("Top 200 by combined score: %d genes\n", nrow(targets_top200)))

# --- Negative control: miRNAs NOT expressed in mESCs -----------------------
# Already-processed, ready to use directly -- same structure as
# targets_all, no parsing needed (unlike the per-family analysis, which is
# pending confirmation of the raw TargetScan .txt format).
targets_negctrl <- readRDS("resources/mouse/Targets_Negative_Control_mESCs.rds")
cat(sprintf("Negative control target list: %d genes\n", nrow(targets_negctrl)))

# --- Our method's activity function, EXACT KO-mean-normalized convention --
# Reproduces separate_KO_control.R's calculate_activity() exactly (KO-group
# mean as the normalization denominator, not the population mean our
# standard Utils.R version uses elsewhere in this project) -- deliberately
# NOT the shared calculate_activity(), since this analysis specifically
# needs the KO-baseline convention.
calculate_activity_ko_normalized <- function(counts_mat, target_ids) {
  cell_sums <- counts_mat[rownames(counts_mat) %in% target_ids, , drop = FALSE] %>% colSums()
  cell_sums_norm <- cell_sums / mean(cell_sums[kos])
  activity <- -log2(cell_sums_norm)
  activity
}

our_all_activity <- calculate_activity_ko_normalized(counts, targets_all$ensembl_gene_id)
our_top200_activity <- calculate_activity_ko_normalized(counts, targets_top200$ensembl_gene_id)

stopifnot("ensembl_gene_id" %in% colnames(targets_negctrl))
our_negctrl_activity <- calculate_activity_ko_normalized(counts, targets_negctrl$ensembl_gene_id)

# --- Per-family target lists (from preprocess_family_targets.R) ------------
targets_mir17 <- readRDS("resources/mouse/Targets_MIR17_mESCs.rds")
targets_mir291 <- readRDS("resources/mouse/Targets_MIR291_mESCs.rds")
targets_mir292 <- readRDS("resources/mouse/Targets_MIR292_mESCs.rds")
cat(sprintf("Per-family target lists: MIR17=%d, MIR291=%d, MIR292=%d genes\n",
            nrow(targets_mir17), nrow(targets_mir291), nrow(targets_mir292)))

our_mir17_activity <- calculate_activity_ko_normalized(counts, targets_mir17$ensembl_gene_id)
our_mir291_activity <- calculate_activity_ko_normalized(counts, targets_mir291$ensembl_gene_id)
our_mir292_activity <- calculate_activity_ko_normalized(counts, targets_mir292$ensembl_gene_id)

# --- Compute, report, plot, and save -- looped over all six variants -------
variants <- list(
  all      = our_all_activity,
  top200   = our_top200_activity,
  negctrl  = our_negctrl_activity,
  MIR17    = our_mir17_activity,
  MIR291   = our_mir291_activity,
  MIR292   = our_mir292_activity
)
variant_labels <- list(
  all      = "Our method, all targets",
  top200   = "Our method, top 200 targets",
  negctrl  = "Our method, NEGATIVE CONTROL (miRNAs not in mESCs)",
  MIR17    = "Our method, MIR-17 family only",
  MIR291   = "Our method, MIR-291 family only",
  MIR292   = "Our method, MIR-292 family only"
)

sep_results <- list()
cat("\n=== Our method: KO vs Control separation, all variants ===\n")
for (v in names(variants)) {
  sep_results[[v]] <- compute_separation(variants[[v]], pheno)
  report_separation(sep_results[[v]], variant_labels[[v]])
}
cat("\nExpectation: NEGATIVE CONTROL should show near-zero separation (overlap close to 100%,\n",
    "mean difference close to 0) -- if it shows strong separation too, that indicates a confound\n",
    "rather than genuine miR-target-specific signal, worth investigating before trusting any of\n",
    "the family-specific results above.\n",
    "Also worth comparing each individual family's separation against the COMBINED (all-targets)\n",
    "result -- if one family dominates the combined signal while the other two contribute little,\n",
    "that's informative about which of the three families is actually driving the KO/Control split.\n")

for (v in names(variants)) {
  p <- plot_separation(sep_results[[v]], variant_labels[[v]])
  print(p)
  ggsave(file.path(OUT_DIR, sprintf("separation_our_method_%s.pdf", v)), p, width = 7, height = 6)
}

saveRDS(sep_results, file.path(OUT_DIR, "res_our_method_separation.rds"))
cat(sprintf("\nSaved to %s\n", OUT_DIR))
