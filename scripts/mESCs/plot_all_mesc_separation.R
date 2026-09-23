library(tidyverse)

BASE_DIR <- "/cfs/klemming/projects/supr/naiss2024-6-235/miRNA_footprint_new"
OUT_DIR <- file.path(BASE_DIR, "analysis/mESCs")

# ---------------------------------------------------------------------------
# Combines separation results across all four methods. Each method has a
# DIFFERENT set of variants available (our method: 9 incl. top200s;
# bayesReact/miReact: 4, per-family + combined, no negctrl since neither
# has a motif for the negative-control miRNAs; miTEA-HiRes: 2, combined +
# negctrl only). Builds one long-format table from whatever's actually
# there, gracefully skipping any method not yet run -- same principle as
# plot_all_comparisons.R for the HEK293 comparisons.
# ---------------------------------------------------------------------------

safe_read <- function(path) {
  if (!file.exists(path)) {
    cat(sprintf("SKIPPED (not found): %s\n", path))
    return(NULL)
  }
  cat(sprintf("Loaded: %s\n", path))
  readRDS(path)
}

# Each result file is a named list of sep_result objects
# ($overlap, $diff, $d, $n_ko, $n_control) -- flatten into a long tibble.
flatten_results <- function(res_list, method_name) {
  if (is.null(res_list)) return(NULL)
  map_dfr(names(res_list), function(variant) {
    r <- res_list[[variant]]
    tibble(method = method_name, variant = variant,
           overlap = r$overlap, diff = r$diff, d = r$d,
           n_ko = r$n_ko, n_control = r$n_control)
  })
}

our_res <- safe_read(file.path(OUT_DIR, "res_our_method_separation.rds"))
bayesreact_res <- safe_read(file.path(OUT_DIR, "res_bayesreact_separation.rds"))
mireact_res <- safe_read(file.path(OUT_DIR, "res_mireact_separation.rds"))
mitea_res <- safe_read(file.path(OUT_DIR, "res_mitea_separation.rds"))

all_results <- bind_rows(
  flatten_results(our_res, "Our method"),
  flatten_results(bayesreact_res, "bayesReact"),
  flatten_results(mireact_res, "miReact"),
  flatten_results(mitea_res, "miTEA-HiRes")
)

if (nrow(all_results) == 0) stop("No result files found -- check OUT_DIR and that at least one method has been run.")

cat("\n=== All loaded results ===\n")
print(all_results %>% select(method, variant, overlap, diff, d), n = Inf)

# Our method uses "all" for the combined/all-targets result; other three
# use "combined" -- unify naming so the cross-method comparison plot below
# can find the right rows from each method consistently.
all_results <- all_results %>%
  mutate(variant_unified = case_when(
    method == "Our method" & variant == "all" ~ "combined",
    TRUE ~ variant
  ))

base_theme <- theme_bw(base_size = 15) +
  theme(panel.grid.minor = element_blank(), axis.text.x = element_text(angle = 20, hjust = 1))

# --- Plot A: combined (all-targets / sum-of-3-families) result, per method -
plot_a_data <- all_results %>% filter(variant_unified == "combined")
p_combined <- ggplot(plot_a_data, aes(x = reorder(method, overlap), y = overlap, fill = method)) +
  geom_col(width = 0.6) +
  geom_text(aes(label = sprintf("%.1f%%", overlap), y = overlap + 2), size = 4) +
  base_theme + theme(legend.position = "none") +
  labs(x = NULL, y = "Overlap % (lower = better separation)",
       title = "KO vs Control separation: combined result, all methods",
       subtitle = "Each method's full/combined target set -- most directly comparable number across methods")
ggsave(file.path(OUT_DIR, "separation_comparison_combined.pdf"), p_combined, width = 7, height = 6, dpi = 600)

# --- Plot B: negative control, where available ------------------------------
plot_b_data <- all_results %>% filter(str_detect(variant, "negctrl"))
if (nrow(plot_b_data) > 0) {
  p_negctrl <- ggplot(plot_b_data, aes(x = method, y = overlap, fill = method)) +
    geom_col(width = 0.6) +
    geom_hline(yintercept = 100, linetype = "dashed", color = "red") +
    geom_text(aes(label = sprintf("%.1f%%", overlap), y = overlap + 2), size = 4) +
    base_theme + theme(legend.position = "none") +
    ylim(0, 110) +
    labs(x = NULL, y = "Overlap %",
         title = "Negative control: expected near 100% (dashed line) if no confound",
         subtitle = "Only methods where a negative-control target/motif was available")
  ggsave(file.path(OUT_DIR, "separation_comparison_negctrl.pdf"), p_negctrl, width = 6, height = 6, dpi = 600)
} else {
  cat("No negative control results found for any method -- skipping that plot.\n")
}

# --- Plot C: per-family breakdown, where available ---------------------------
plot_c_data <- all_results %>% filter(variant_unified %in% c("MIR17", "MIR291", "MIR292"))
if (nrow(plot_c_data) > 0) {
  p_family <- ggplot(plot_c_data, aes(x = variant_unified, y = overlap, fill = method)) +
    geom_col(position = position_dodge(width = 0.7), width = 0.6) +
    geom_text(aes(label = sprintf("%.0f%%", overlap)), position = position_dodge(width = 0.7),
               vjust = -0.3, size = 3) +
    base_theme +
    labs(x = NULL, y = "Overlap % (lower = better separation)", fill = NULL,
         title = "Per-family separation, all methods with per-family results available")
  ggsave(file.path(OUT_DIR, "separation_comparison_per_family.pdf"), p_family, width = 8, height = 6, dpi = 600)
} else {
  cat("No per-family results found for any method -- skipping that plot.\n")
}

print(p_combined)
if (exists("p_negctrl")) print(p_negctrl)
if (exists("p_family")) print(p_family)

saveRDS(all_results, file.path(OUT_DIR, "res_all_methods_separation_combined.rds"))
cat(sprintf("\nSaved combined table and plots to %s\n", OUT_DIR))
