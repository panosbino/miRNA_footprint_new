library(tidyverse)

BASE_DIR <- "/cfs/klemming/projects/supr/naiss2024-6-235/miRNA_footprint_new"
OUT_DIR <- file.path(BASE_DIR, "analysis/mESCs")
MITEA_INPUT_DIR <- file.path(OUT_DIR, "mitea_input")

source(file.path(BASE_DIR, "scripts/mESCs/mESC_separation_utils.R"))

setwd(BASE_DIR)
pheno <- read.delim("./datasets/mESCs/phenotype.csv", sep = ",")
pheno$Exp[pheno$Exp == "WT"] <- "Control"

load_as_named_vec <- function(csv_file) {
  df <- read.csv(file.path(MITEA_INPUT_DIR, csv_file))
  setNames(df$activity, df$Cell_ID)
}

activity_combined <- load_as_named_vec("mitea_scores_combined.csv")
activity_negctrl <- load_as_named_vec("mitea_scores_negctrl.csv")

variants <- list(combined = activity_combined, negctrl = activity_negctrl)
variant_labels <- list(combined = "miTEA-HiRes, combined (3-family target list)",
                        negctrl = "miTEA-HiRes, NEGATIVE CONTROL")

sep_results <- list()
cat("\n=== miTEA-HiRes: KO vs Control separation ===\n")
for (v in names(variants)) {
  sep_results[[v]] <- compute_separation(variants[[v]], pheno)
  report_separation(sep_results[[v]], variant_labels[[v]])
}
cat("\nExpectation: negative control should show near-zero separation -- check this before trusting the combined result above.\n")

for (v in names(variants)) {
  p <- plot_separation(sep_results[[v]], variant_labels[[v]])
  print(p)
  ggsave(file.path(OUT_DIR, sprintf("separation_mitea_%s.pdf", v)), p, width = 7, height = 6)
}

saveRDS(sep_results, file.path(OUT_DIR, "res_mitea_separation.rds"))
cat(sprintf("\nSaved to %s\n", OUT_DIR))
