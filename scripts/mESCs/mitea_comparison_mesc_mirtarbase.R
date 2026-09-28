library(tidyverse)

BASE_DIR <- "/cfs/klemming/projects/supr/naiss2024-6-235/miRNA_footprint_new"
OUT_DIR <- file.path(BASE_DIR, "analysis/mESCs")
MITEA_INPUT_DIR <- file.path(OUT_DIR, "mitea_input")

source(file.path(BASE_DIR, "scripts/mESCs/mESC_separation_utils.R"))

setwd(BASE_DIR)
pheno <- read.delim("./datasets/mESCs/phenotype.csv", sep = ",")
pheno$Exp[pheno$Exp == "WT"] <- "Control"

# Target counts written by run_mitea_mesc_mirtarbase.py -- reported next to
# every result because validated mouse targets for these miRNAs are few
# (measured from the bundled file: ~48 genes MIR17, ~11 MIR291, 2 MIR292),
# which directly limits how far each result can be trusted.
counts_path <- file.path(MITEA_INPUT_DIR, "mirtarbase_target_counts.csv")
target_counts <- if (file.exists(counts_path)) read.csv(counts_path) else NULL
if (is.null(target_counts)) cat("NOTE: mirtarbase_target_counts.csv not found; target counts will not be shown.\n")

load_as_named_vec <- function(variant) {
  f <- file.path(MITEA_INPUT_DIR, sprintf("mitea_scores_mirtarbase_%s.csv", variant))
  if (!file.exists(f)) {
    cat(sprintf("SKIPPED %s (no score file -- the Python step skips a set with no targets in the data)\n", variant))
    return(NULL)
  }
  df <- read.csv(f)
  setNames(df$activity, df$Cell_ID)
}

variant_names <- c("combined", "MIR17", "MIR291", "MIR292")
variants <- setNames(lapply(variant_names, load_as_named_vec), variant_names)
variants <- variants[!sapply(variants, is.null)]
if (length(variants) == 0) stop("No miRTarBase-based score files found. Run run_mitea_mesc_mirtarbase.py first.")

variant_labels <- list(
  combined = "miTEA-HiRes + miRTarBase, combined (3 families)",
  MIR17    = "miTEA-HiRes + miRTarBase, MIR-17 family",
  MIR291   = "miTEA-HiRes + miRTarBase, MIR-291 family",
  MIR292   = "miTEA-HiRes + miRTarBase, MIR-292 family"
)

sep_results <- list()
cat("\n=== miTEA-HiRes (native mouse miRTarBase): KO vs Control separation ===\n")
for (v in names(variants)) {
  sep_results[[v]] <- compute_separation(variants[[v]], pheno)
  report_separation(sep_results[[v]], variant_labels[[v]])
  if (!is.null(target_counts)) {
    tc <- target_counts[target_counts$variant == v, ]
    if (nrow(tc) == 1) {
      cat(sprintf("    targets: %d unique validated, %d present in expression matrix%s\n",
                  tc$n_unique_targets, tc$n_targets_in_data,
                  ifelse(tc$n_targets_in_data < 10, "   <-- FEW TARGETS: treat as unreliable", "")))
    }
  }
}
cat("\nReading guide: lower overlap = better separation. With so few validated targets, MIR-291 and\n",
    "especially MIR-292 differences from the other variants may reflect noise rather than the tool.\n",
    "Compare against the TargetScan-list miTEA-HiRes run (res_mitea_separation.rds) and the other\n",
    "methods with the same caveat in mind.\n")

for (v in names(variants)) {
  p <- plot_separation(sep_results[[v]], variant_labels[[v]])
  print(p)
  ggsave(file.path(OUT_DIR, sprintf("separation_mitea_mirtarbase_%s.pdf", v)), p, width = 7, height = 6)
}

saveRDS(sep_results, file.path(OUT_DIR, "res_mitea_mirtarbase_separation.rds"))
cat(sprintf("\nSaved to %s\n", OUT_DIR))
