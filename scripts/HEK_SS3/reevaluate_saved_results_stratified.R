library(tidyverse)

# ==============================================================================
# Re-score every tool's SAVED per-cell results against
#   {GFP mRNA, GFP fluorescence} x {all cells, Dox-induced cells only}
# WITHOUT rerunning miReact / bayesReact / miTEA-HiRes.
#
# Why this works: each res_*_comparison_scran.rds already stores per-cell scores
# for exactly the QC-passing cells (GFP mRNA < 99th pct), which is the same cell
# set load_gfp_truth() uses. Only the ground truth and the cell subsetting change.
#
# Run the per-tool scripts instead only if the tool output itself must change.
# ==============================================================================

BASE_DIR      <- Sys.getenv("MIRNA_BASE_DIR", "/cfs/klemming/projects/supr/naiss2024-6-235/miRNA_footprint_new")
PROCESSED_DIR <- file.path(BASE_DIR, "datasets/HEK_SS3/processed")
OUT_DIR       <- Sys.getenv("MIRNA_COMPARISON_DIR", file.path(BASE_DIR, "analysis/HEK_SS3/comparisons"))
GFP_UPPER_PCT <- 0.99
N_BOOT        <- 1000
OURS_LABEL    <- "Our method (TargetScan top200)"

source(file.path(BASE_DIR, "scripts/Utils.R"))
gfp_truth <- load_gfp_truth(PROCESSED_DIR, upper_pct = GFP_UPPER_PCT)

# file -> list(element in the rds = method label). Our method (TargetScan top200)
# is stored in every file; we keep one copy and check the others agree.
SOURCES <- list(
  res_mireact_comparison_scran.rds = c(our_targetscan = OURS_LABEL,
                                       our_tarbase    = "Our method (TarBase)",
                                       mireact        = "miReact (TarBase)"),
  res_bayesreact_comparison_scran.rds = c(our = OURS_LABEL, bayesreact = "bayesReact (motif)"),
  res_mireact_default_comparison_scran.rds = c(our = OURS_LABEL, mireact_default = "miReact (motif)"),
  res_mitea_comparison_scran.rds = c(our = OURS_LABEL,
                                     mitea_mirtarbase        = "miTEA-HiRes (miRTarBase)",
                                     mitea_tarbase           = "miTEA-HiRes (TarBase)",
                                     mitea_targetscan_all    = "miTEA-HiRes (TargetScan (all))",
                                     mitea_targetscan_top200 = "miTEA-HiRes (TargetScan (top 200))")
)

# The per-cell data frames hold Cell_ID, one score column, and the old GFP column (eGFP).
extract_scores <- function(df, where) {
  if (!is.data.frame(df) || !"Cell_ID" %in% colnames(df)) stop("No per-cell data frame with Cell_ID at ", where)
  candidates <- setdiff(colnames(df)[vapply(df, is.numeric, logical(1))], c("eGFP", "GFP", "gfp_mrna", "gfp_facs"))
  if (length(candidates) != 1)
    stop(sprintf("Can't identify the score column at %s (numeric columns: %s)", where, paste(candidates, collapse = ", ")))
  tibble(Cell_ID = df$Cell_ID, score = df[[candidates]])
}

scores <- list()
for (f in names(SOURCES)) {
  path <- file.path(OUT_DIR, f)
  if (!file.exists(path)) { cat(sprintf("SKIPPED (not found): %s\n", path)); next }
  obj <- readRDS(path)
  cat(sprintf("Loaded: %s\n", f))
  for (el in names(SOURCES[[f]])) {
    label <- SOURCES[[f]][[el]]
    if (is.null(obj[[el]])) { cat(sprintf("  WARNING: '%s' not in %s -- skipped\n", el, f)); next }
    s <- extract_scores(obj[[el]], paste0(f, "$", el))
    if (label %in% names(scores)) {
      # Consistency check: the same method stored in two files must agree.
      m <- inner_join(scores[[label]], s, by = "Cell_ID")
      max_diff <- max(abs(m$score.x - m$score.y))
      cat(sprintf("  consistency: %s in %s vs first copy: %d shared cells, max |diff| = %.2e%s\n",
                  label, f, nrow(m), max_diff, if (max_diff > 1e-8) "   *** MISMATCH -- runs used different inputs ***" else ""))
    } else {
      scores[[label]] <- s
    }
  }
}
if (length(scores) == 0) stop("No saved results found in ", OUT_DIR)

cat("\n=== Stratified evaluation ===\n")
stratified <- imap_dfr(scores, function(s, label)
  evaluate_score_vs_gfp(s, "score", label, gfp_truth, n_boot = N_BOOT))

cat("\n=== Paired differences vs our method (same cells) ===\n")
paired <- if (OURS_LABEL %in% names(scores)) {
  imap_dfr(scores[names(scores) != OURS_LABEL], function(s, label)
    paired_rho_difference(scores[[OURS_LABEL]], "score", OURS_LABEL, s, "score", label,
                          gfp_truth, n_boot = N_BOOT))
} else tibble()

write.csv(stratified, file.path(OUT_DIR, "stratified_all_tools.csv"), row.names = FALSE)
write.csv(paired, file.path(OUT_DIR, "paired_vs_ours_all_tools.csv"), row.names = FALSE)
p <- plot_stratified(stratified)
ggsave(file.path(OUT_DIR, "stratified_all_tools.pdf"), p, width = 12, height = 1.5 + 0.9 * length(scores))
ggsave(file.path(OUT_DIR, "stratified_all_tools.png"), p, width = 12, height = 1.5 + 0.9 * length(scores), dpi = 200)
cat(sprintf("\nWrote stratified_all_tools.{csv,pdf,png} and paired_vs_ours_all_tools.csv to %s\n", OUT_DIR))
