library(tidyverse)

BASE_DIR <- "/cfs/klemming/projects/supr/naiss2024-6-235/miRNA_footprint_new"
OUT_DIR <- file.path(BASE_DIR, "analysis/mESCs")

source(file.path(BASE_DIR, "scripts/mESCs/mESC_separation_utils.R"))

setwd(BASE_DIR)
pheno <- read.delim("./datasets/mESCs/phenotype.csv", sep = ",")
pheno$Exp[pheno$Exp == "WT"] <- "Control"

tracking_path <- file.path(OUT_DIR, "mireact_mesc_job_tracking.rds")
if (!file.exists(tracking_path)) stop(sprintf("No tracking file found at %s. Run launch_mireact_mesc.R first.", tracking_path))
tracking_info <- readRDS(tracking_path)
cat(sprintf("Job submitted at: %s\nExpected directory: %s\n", tracking_info$submitted_at, tracking_info$wd))

result_file <- file.path(tracking_info$wd, tracking_info$out_file)
if (!file.exists(result_file)) {
  stop(sprintf("Result file not found at %s -- job may still be running (check squeue -u $USER) or may have failed (check %s/Rscript-*.out).",
               result_file, tracking_info$wd))
}
cat(sprintf("Found result file: %s\n", result_file))

ma <- readRDS(result_file)
cat(sprintf("Motif-activity matrix: %d motifs x %d cells\n", nrow(ma), ncol(ma)))
cat("Sample column names:\n"); print(head(colnames(ma)))

MOTIFS <- list(MIR17 = "GCACTTT", MIR291 = "AGCACTT", MIR292 = "GGCACTT")

# Mandatory sanity check: confirm each motif is actually a valid row key --
# same pattern as the HEK293 default-mode collection script.
for (fam in names(MOTIFS)) {
  if (!(MOTIFS[[fam]] %in% rownames(ma))) {
    stop(sprintf("'%s' (%s) not found in rownames(ma). Check rownames(ma) directly to see the actual format.",
                 MOTIFS[[fam]], fam))
  }
}
cat("*** Confirmed all three motifs present as rows in ma. ***\n")

activity_mir17 <- setNames(ma[MOTIFS$MIR17, ], colnames(ma))
activity_mir291 <- setNames(ma[MOTIFS$MIR291, ], colnames(ma))
activity_mir292 <- setNames(ma[MOTIFS$MIR292, ], colnames(ma))
activity_combined <- activity_mir17 + activity_mir291 + activity_mir292

variants <- list(MIR17 = activity_mir17, MIR291 = activity_mir291,
                  MIR292 = activity_mir292, combined = activity_combined)
variant_labels <- list(MIR17 = "miReact, MIR-17 family (motif mode)",
                        MIR291 = "miReact, MIR-291 family (motif mode)",
                        MIR292 = "miReact, MIR-292 family (motif mode)",
                        combined = "miReact, combined (sum of 3 families, motif mode)")

sep_results <- list()
cat("\n=== miReact: KO vs Control separation ===\n")
for (v in names(variants)) {
  sep_results[[v]] <- compute_separation(variants[[v]], pheno)
  report_separation(sep_results[[v]], variant_labels[[v]])
}

for (v in names(variants)) {
  p <- plot_separation(sep_results[[v]], variant_labels[[v]])
  print(p)
  ggsave(file.path(OUT_DIR, sprintf("separation_mireact_%s.pdf", v)), p, width = 7, height = 6)
}

saveRDS(sep_results, file.path(OUT_DIR, "res_mireact_separation.rds"))
cat(sprintf("\nSaved to %s\n", OUT_DIR))
