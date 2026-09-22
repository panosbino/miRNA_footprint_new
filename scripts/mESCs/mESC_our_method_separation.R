library(tidyverse)

BASE_DIR <- "/cfs/klemming/projects/supr/naiss2024-6-235/miRNA_footprint_new"
OUT_DIR <- file.path(BASE_DIR, "analysis/mESCs/KO_control_separation")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

source(file.path(BASE_DIR, "scripts/HEK_SS3/comparisons/mESC_separation_utils.R"))

# --- Data (exactly matching separate_KO_control.R's loading) ---------------
setwd(BASE_DIR)
pheno <- read.delim("./datasets/mESCs/phenotype.csv", sep = ",")
pheno$Exp[pheno$Exp == "WT"] <- "Control"
kos <- pheno[pheno$Exp == "KO", ]$cells

counts <- readRDS("datasets/mESCs/final_counts.rds")
counts <- counts[, colnames(counts) %in% pheno$cells]

targets_all <- readRDS("resources/Targets_combined_top3_families_mESCs.rds")
targets_all <- targets_all %>%
  group_by(ensembl_gene_id) %>%
  mutate(score = sum(Cumulative.weighted.context...score)) %>%
  ungroup() %>%
  distinct(pick(ensembl_gene_id), .keep_all = TRUE)

cat(sprintf("Combined top-3-families target list: %d unique genes\n", nrow(targets_all)))

targets_top200 <- targets_all %>% arrange(score) %>% slice_head(n = 200)  # most negative/strongest score first, matching TargetScan convention used throughout this project
cat(sprintf("Top 200 by combined score: %d genes\n", nrow(targets_top200)))

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

sep_all <- compute_separation(our_all_activity, pheno)
sep_top200 <- compute_separation(our_top200_activity, pheno)

cat("\n=== Our method: KO vs Control separation ===\n")
report_separation(sep_all, "Our method, all targets")
report_separation(sep_top200, "Our method, top 200 targets")

p_all <- plot_separation(sep_all, "Our method (all targets)")
p_top200 <- plot_separation(sep_top200, "Our method (top 200 targets)")
print(p_all)
print(p_top200)

ggsave(file.path(OUT_DIR, "separation_our_method_all.pdf"), p_all, width = 7, height = 6)
ggsave(file.path(OUT_DIR, "separation_our_method_top200.pdf"), p_top200, width = 7, height = 6)

saveRDS(list(all = sep_all, top200 = sep_top200,
             targets_all = targets_all$ensembl_gene_id, targets_top200 = targets_top200$ensembl_gene_id),
        file.path(OUT_DIR, "res_our_method_separation.rds"))

cat(sprintf("\nSaved to %s\n", OUT_DIR))
