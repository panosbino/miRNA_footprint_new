library(dplyr)
library(stringr)
library(biomaRt)

# ---------------------------------------------------------------------------
# RUN THIS ON THE LOGIN NODE, NOT VIA SBATCH. Uses a live biomaRt query
# (transcript ID -> gene ID mapping) -- compute nodes on Dardel typically
# lack outbound internet access, the same issue we already hit and worked
# around for the human protein-coding gene list earlier in this project.
# This is a ONE-TIME step; the resulting RDS files are then just loaded
# directly by the main analysis script, same "cache the network-dependent
# step separately" pattern used there.
#
# Reuses the EXACT transformation logic from preprocess_targets.R
# (strip transcript version -> biomaRt transcript-to-gene mapping -> sort
# by score -> drop unmapped genes), applied to each of the three family
# files SEPARATELY (no common_cols intersection needed here, since we're
# not combining files the way the original combined-target preprocessing
# did).
# ---------------------------------------------------------------------------

BASE_DIR <- "/cfs/klemming/projects/supr/naiss2024-6-235/miRNA_footprint_new"
setwd(BASE_DIR)

# useEnsembl(), not the deprecated useMart() -- same fix already established
# earlier in this project for exactly this class of connection issue.
ensembl_mm <- useEnsembl(biomart = "ensembl", dataset = "mmusculus_gene_ensembl")

first.word <- function(my.string) {
  unlist(str_split(my.string, fixed(".")))[1]
}

switch_names_targets <- function(targets) {
  no_genes <- dim(targets)[1]
  genes <- getBM(attributes = c('ensembl_gene_id', 'ensembl_transcript_id'),
                  filters = 'ensembl_transcript_id',
                  values = targets$Representative.transcript %>% unique(),
                  mart = ensembl_mm)
  targets2 <- merge(targets, genes, by.x = "Representative.transcript", by.y = "ensembl_transcript_id", all.x = TRUE)
  cat(dim(genes)[1], "ENSEMBL transcript IDs successfully changed ENSEMBL gene IDs, out of", no_genes, "total genes\n")
  return(targets2)
}

preprocess_one_family <- function(txt_path) {
  targets <- read.delim(txt_path, sep = "\t")
  targets$Representative.transcript <- sapply(targets$Representative.transcript, first.word)
  targets_new <- switch_names_targets(targets = targets)
  targets_new <- arrange(targets_new, Cumulative.weighted.context...score)
  targets_new <- targets_new[!is.na(targets_new$ensembl_gene_id), ]
  targets_new
}

families <- list(
  MIR17  = "resources/mouse/Targetscan_files_mouse/mouse/TargetScan8.0__miR-17-5p_20-5p_93-5p_106-5p.predicted_targets.txt",
  MIR291 = "resources/mouse/Targetscan_files_mouse/mouse/TargetScan8.0__miR-291-3p_294-3p_295-3p_302-3p.predicted_targets.txt",
  MIR292 = "resources/mouse/Targetscan_files_mouse/mouse/TargetScan8.0__miR-292a-3p_467a-5p.predicted_targets.txt"
)

for (fam_name in names(families)) {
  cat(sprintf("\n=== Processing %s ===\n", fam_name))
  targets_processed <- preprocess_one_family(families[[fam_name]])
  out_path <- sprintf("resources/mouse/Targets_%s_mESCs.rds", fam_name)
  saveRDS(targets_processed, out_path)
  cat(sprintf("Wrote %s (%d genes)\n", out_path, nrow(targets_processed)))
}

cat("\nDone. Three files written to resources/mouse/: Targets_MIR17_mESCs.rds, Targets_MIR291_mESCs.rds, Targets_MIR292_mESCs.rds\n")
