library(tidyverse)

#BASE_DIR <- "/cfs/klemming/projects/supr/naiss2024-6-235/miRNA_footprint_new"
BASE_DIR <- Sys.getenv("MIRNA_BASE_DIR", "~/Desktop/Projects/miRNA_footprint_new/")
source(file.path(BASE_DIR, "scripts/Utils.R"))   # for plot_stratified()

OUT_DIR <- file.path(BASE_DIR, "analysis/HEK_SS3/comparisons")

# ---------------------------------------------------------------------------
# Robust to not knowing the exact internal field names in each RDS (this
# conversation has spanned a long time and field names may have drifted
# slightly from what was originally documented). Rather than hard-code
# names that might be wrong, auto-detect every cor.test-like object inside
# each loaded list -- anything with BOTH $estimate and $p.value is treated
# as a correlation result, labeled by its list name.
# ---------------------------------------------------------------------------

safe_read <- function(path) {
  if (!file.exists(path)) {
    cat(sprintf("SKIPPED (not found): %s\n", path))
    return(NULL)
  }
  cat(sprintf("Loaded: %s\n", path))
  readRDS(path)
}

is_cor_test_like <- function(x) {
  # Skip data frames (e.g. the new stratified/paired tibbles): probing them with $estimate warns.
  is.list(x) && !is.data.frame(x) && !is.null(x$estimate) && !is.null(x$p.value) && is.numeric(x$estimate)
}

extract_all_cors <- function(rds_obj, source_label) {
  if (is.null(rds_obj)) return(NULL)
  cor_fields <- names(rds_obj)[sapply(rds_obj, is_cor_test_like)]
  if (length(cor_fields) == 0) {
    cat(sprintf("  WARNING: no cor.test-like objects auto-detected in %s -- check its structure manually with names(readRDS(...)).\n", source_label))
    return(NULL)
  }
  cat(sprintf("  Found %d correlation result(s) in %s: %s\n", length(cor_fields), source_label, paste(cor_fields, collapse = ", ")))
  map_dfr(cor_fields, function(fn) {
    ct <- rds_obj[[fn]]
    df_name <- sub("_cor$", "", fn)
    n_val <- if (!is.null(rds_obj[[df_name]]) && is.data.frame(rds_obj[[df_name]])) nrow(rds_obj[[df_name]]) else NA
    tibble(source = source_label, result_type = fn, rho = unname(ct$estimate), p_value = ct$p.value, n = n_val)
  })
}

mireact_tarbase <- safe_read(file.path(OUT_DIR, "res_mireact_comparison_scran.rds"))
bayesreact_res <- safe_read(file.path(OUT_DIR, "res_bayesreact_comparison_scran.rds"))
mireact_default <- safe_read(file.path(OUT_DIR, "res_mireact_default_comparison_scran.rds"))
mitea_res <- safe_read(file.path(OUT_DIR, "res_mitea_comparison_scran.rds"))

cat("\n=== Auto-detecting correlation results in each file ===\n")
all_results <- bind_rows(
  extract_all_cors(mireact_tarbase, "mireact_comparison (TarBase mode)"),
  extract_all_cors(bayesreact_res, "bayesreact_comparison"),
  extract_all_cors(mireact_default, "mireact_default_comparison (motif mode)"),
  extract_all_cors(mitea_res, "mitea_comparison")
)

if (nrow(all_results) == 0) stop("No correlation results found in any file -- check OUT_DIR and file contents directly.")

cat("\n=== All extracted results ===\n")
print(all_results, n = Inf)

label_map <- c(
  our_targetscan_cor = "Our method\n(TargetScan)",
  our_tarbase_cor = "Our method\n(TarBase)",
  our_cor = "Our method\n(TargetScan)",
  mireact_cor = "miReact\n(TarBase)",
  mireact_default_cor = "miReact\n(motif mode)",
  bayesreact_cor = "bayesReact\n(motif mode)",
  mitea_cor_mirtarbase = "miTEA-HiRes\n(miRTarBase)",
  mitea_cor_tarbase = "miTEA-HiRes\n(TarBase)",
  mitea_cor_targetscan_all = "miTEA-HiRes\n(TargetScan, all)",
  mitea_cor_targetscan_top200 = "miTEA-HiRes\n(TargetScan, top200)"
)

all_results <- all_results %>%
  mutate(label = ifelse(result_type %in% names(label_map), label_map[result_type], result_type))

y_limits <- c(min(0, floor(min(all_results$rho, na.rm = TRUE) * 10) / 10),
              max(1, ceiling(max(all_results$rho, na.rm = TRUE) * 10) / 10))

# "Our method" is stored in every file under the same label; without de-duplication
# geom_col stacks the copies (~3.07 tall) and ylim() silently drops the excess.
plot_data <- all_results %>% distinct(label, .keep_all = TRUE)
p_all <- ggplot(plot_data, aes(x = reorder(label, rho), y = rho, fill = source)) +
  geom_col(width = 0.65) +
  geom_text(aes(label = sprintf("%.3f", rho), y = rho + 0.03), size = 3.8) +
  ylim(y_limits[1], y_limits[2]) +
  theme_bw(base_size = 14) +
  theme(panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 30, hjust = 1),
        legend.position = "bottom", legend.title = element_blank()) +
  labs(x = NULL, y = "Spearman rho",
       title = "Tool comparison, HEK293 dataset (scran-normalized)",
       subtitle = "Every correlation result found across all four comparison scripts")

print(p_all)
ggsave(file.path(OUT_DIR, "comparison_all_results.pdf"), p_all, width = 10, height = 7, dpi = 600)

default_types <- c("our_targetscan_cor", "our_cor", "mireact_cor", "mireact_default_cor",
                    "bayesreact_cor", "mitea_cor_mirtarbase")
default_data <- all_results %>% filter(result_type %in% default_types) %>% distinct(label, .keep_all = TRUE)

if (nrow(default_data) > 0) {
  p_default <- ggplot(default_data, aes(x = reorder(label, rho), y = rho, fill = label)) +
    geom_col(width = 0.6) +
    geom_text(aes(label = sprintf("%.3f", rho), y = rho + 0.03), size = 4) +
    ylim(y_limits[1], y_limits[2]) +
    theme_bw(base_size = 15) +
    theme(panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 20, hjust = 1),
          legend.position = "none") +
    labs(x = NULL, y = "Spearman rho",
         title = "Each method, as normally used (HEK293, scran-normalized)")
  print(p_default)
  ggsave(file.path(OUT_DIR, "comparison_default_per_method.pdf"), p_default, width = 8, height = 6, dpi = 600)
}

saveRDS(all_results, file.path(OUT_DIR, "res_all_tools_comparison_combined.rds"))
cat(sprintf("\nSaved plots and combined table to %s\n", OUT_DIR))

# ---------------------------------------------------------------------------
# Stratified results ({GFP mRNA, fluorescence} x {all, induced}), present only
# in .rds files written by the updated tool scripts. For results from older
# runs, use reevaluate_saved_results_stratified.R instead (no tool reruns).
# ---------------------------------------------------------------------------
stratified_all <- bind_rows(lapply(list(mireact_tarbase, bayesreact_res, mireact_default, mitea_res),
                                   function(x) if (!is.null(x)) x$stratified))
if (nrow(stratified_all) > 0) {
  # "Our method" is evaluated in every script on the same cells -> keep one copy.
  stratified_all <- distinct(stratified_all, method, readout, cell_set, .keep_all = TRUE)
  p_strat <- plot_stratified(stratified_all)
  print(p_strat)
  n_methods <- n_distinct(stratified_all$method)
  ggsave(file.path(OUT_DIR, "comparison_stratified.pdf"), p_strat, width = 12, height = 1.5 + 0.9 * n_methods)
  write.csv(stratified_all, file.path(OUT_DIR, "comparison_stratified.csv"), row.names = FALSE)
  cat(sprintf("Saved stratified comparison (%d methods) to %s\n", n_methods, OUT_DIR))
} else {
  cat("No stratified results inside the .rds files (older runs) -- run reevaluate_saved_results_stratified.R.\n")
}
