library(tidyverse)

# ---------------------------------------------------------------------------
# bayesReact, mESC Control vs KO. Reuses every established fix from the
# HEK293 rebuild (process_raw_input 50%-threshold bypass, build_seq_list
# dedup, single-motif output structure) -- only genuinely new pieces here:
# mouse sequence data, three separate motif runs summed into one score
# (per explicit instruction), and the KO-mean separation-statistic
# framework instead of GFP correlation.
# ---------------------------------------------------------------------------

BASE_DIR <- "/cfs/klemming/projects/supr/naiss2024-6-235/miRNA_footprint_new"
TOOLS_DIR <- file.path(BASE_DIR, "tools")
MIREACT_DIR <- file.path(TOOLS_DIR, "miReact")
OUT_DIR <- file.path(BASE_DIR, "analysis/mESCs")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

.libPaths(c(file.path(TOOLS_DIR, "R_library"), .libPaths()))
library(bayesReact)

source(file.path(BASE_DIR, "scripts/mESCs/mESC_separation_utils.R"))

# Three target-site motifs (already reverse-complemented from the raw miRNA
# seeds -- see earlier verification), one per family, per explicit
# confirmation these are correct.
MOTIFS <- list(MIR17 = "GCACTTT", MIR291 = "AGCACTT", MIR292 = "GGCACTT")

setwd(BASE_DIR)
pheno <- read.delim("./datasets/mESCs/phenotype.csv", sep = ",")
pheno$Exp[pheno$Exp == "WT"] <- "Control"
kos <- pheno[pheno$Exp == "KO", ]$cells

counts <- readRDS("datasets/mESCs/final_counts.rds")
counts <- as.matrix(counts[, colnames(counts) %in% pheno$cells])
cat(sprintf("Counts: %d genes x %d cells\n", nrow(counts), ncol(counts)))

# --- Mouse sequence data, deduplicated (same fix as the HEK293 build) ------
raw_seqs <- readRDS(file.path(MIREACT_DIR, "seqs", "mm.utr3.seqs.rds"))
seqs_deduped <- raw_seqs %>%
  group_by(gid) %>%
  filter(nchar == max(nchar)) %>%
  slice(1) %>%
  ungroup() %>%
  select(gid, sequence, nchar) %>%
  as.data.frame()
stopifnot(anyDuplicated(seqs_deduped$gid) == 0)
cat(sprintf("Deduplicated mouse sequence data: %d genes\n", nrow(seqs_deduped)))

# --- Shared prep (identical for every motif, so run once) ------------------
dir.create(file.path(OUT_DIR, "bayesReact_out"), recursive = TRUE, showWarnings = FALSE)
exp_norm <- bayesReact::norm_scale_seq(counts, data_type = "CPM", save_rds = FALSE)
gene_set <- intersect(seqs_deduped$gid, rownames(exp_norm))
cat(sprintf("Matched genes: %d\n", length(gene_set)))
exp_matched <- exp_norm[gene_set, , drop = FALSE]
seqs_matched <- seqs_deduped[seqs_deduped$gid %in% gene_set, ]
seqs_matched <- seqs_matched[match(gene_set, seqs_matched$gid), ]
seqlist_out <- bayesReact::build_seq_list(seqs_matched, gene_id = "gid")
FC_rank_path <- bayesReact::rank_seq(exp_matched, data_type = "norm_scale_exp",
                                      path = file.path(OUT_DIR, "bayesReact_out/"))

# --- Run bayesReact once per family motif -----------------------------------
run_one_motif <- function(motif_seq, fam_name) {
  cat(sprintf("\n=== bayesReact: %s (motif %s) ===\n", fam_name, motif_seq))
  motif_paths <- bayesReact::motif_prob(motif_seq, seqlist_out$seqs, seqlist_out$seqlist,
                                         paths = FALSE, cores = parallel::detectCores(),
                                         out_path = file.path(OUT_DIR, "bayesReact_out/"), include_counts = TRUE)
  result <- bayesReact_core(
    lst_data = list(FC_rank = FC_rank_path,
                     motif_probs = motif_paths$motif_probs_path,
                     motif_counts = motif_paths$motif_counts_path),
    model = "bayesReact", output_type = "activity_summary", CI = c(0.1, 0.9)
  )
  activity <- result$activity
  names(activity) <- rownames(result)
  cat(sprintf("MCMC diagnostics for %s: median n_eff = %.0f, max Rhat = %.4f\n",
              fam_name, median(result$n_eff, na.rm = TRUE), max(result$Rhat, na.rm = TRUE)))
  activity
}

activity_mir17 <- run_one_motif(MOTIFS$MIR17, "MIR17")
activity_mir291 <- run_one_motif(MOTIFS$MIR291, "MIR291")
activity_mir292 <- run_one_motif(MOTIFS$MIR292, "MIR292")

# --- Sum the three, per explicit instruction --------------------------------
stopifnot(identical(names(activity_mir17), names(activity_mir291)),
          identical(names(activity_mir17), names(activity_mir292)))
activity_combined <- activity_mir17 + activity_mir291 + activity_mir292

# --- Separation analysis: each family individually + combined --------------
variants <- list(MIR17 = activity_mir17, MIR291 = activity_mir291,
                  MIR292 = activity_mir292, combined = activity_combined)
variant_labels <- list(MIR17 = "bayesReact, MIR-17 family",
                        MIR291 = "bayesReact, MIR-291 family",
                        MIR292 = "bayesReact, MIR-292 family",
                        combined = "bayesReact, combined (sum of 3 families)")

sep_results <- list()
cat("\n=== bayesReact: KO vs Control separation ===\n")
for (v in names(variants)) {
  sep_results[[v]] <- compute_separation(variants[[v]], pheno)
  report_separation(sep_results[[v]], variant_labels[[v]])
}

for (v in names(variants)) {
  p <- plot_separation(sep_results[[v]], variant_labels[[v]])
  print(p)
  ggsave(file.path(OUT_DIR, sprintf("separation_bayesreact_%s.pdf", v)), p, width = 7, height = 6)
}

saveRDS(sep_results, file.path(OUT_DIR, "res_bayesreact_separation.rds"))
cat(sprintf("\nSaved to %s\n", OUT_DIR))
