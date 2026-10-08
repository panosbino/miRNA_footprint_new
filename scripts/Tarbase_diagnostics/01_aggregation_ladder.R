# ==============================================================================
# 01  Aggregation ladder
#
# Changes our method one step at a time toward a rank-based method, on both the
# TarBase and TargetScan top-200 lists:
#   A raw sum -> B equal weight -> C log ratio -> C2 background-corrected -> D rank
#
# Reading the result: the step where TarBase performance jumps identifies the
# cause. Jump at B = weighting by expression; at C = outliers / scale; at C2 =
# a cell-wide shift shared by all genes; at D = rank robustness beyond C.
# If TargetScan stays flat across steps while TarBase climbs, the problem is
# specific to long, unranked lists.
#
# Each list x variant also gets a size-matched random-gene null (median and 95th
# percentile over DIAG_N_RANDOM draws), so a rise can be checked against what
# any gene set of that size would give.
#
# Outputs (analysis/HEK_SS3/tarbase_diagnostics/):
#   01_ladder_stratified.csv   rho + 95% CI per list x variant x readout x cell set
#   01_ladder_null.csv         random-gene null per list size x variant
#   01_ladder.pdf              figure
# ==============================================================================

source(file.path(Sys.getenv("MIRNA_BASE_DIR", getwd()), "scripts/HEK_SS3/tarbase_diagnostics/tarbase_diag_common.R"))

cat("=== Ladder: stratified evaluation with bootstrap CIs ===\n")
ladder <- map_dfr(names(LISTS), function(list_name) {
  map_dfr(names(VARIANTS), function(v) {
    evaluate(score_variant(LISTS[[list_name]], v), sprintf("%s | %s", list_name, v)) |>
      mutate(list = list_name, variant = v, variant_label = VARIANTS[[v]], n_genes = length(LISTS[[list_name]]))
  })
})
write.csv(ladder, file.path(OUT_DIR, "01_ladder_stratified.csv"), row.names = FALSE)

cat(sprintf("\n=== Random-gene null: %d draws per list size x variant (induced cells, fluorescence) ===\n", N_RANDOM))
null <- map_dfr(names(LISTS), function(list_name) {
  n <- length(LISTS[[list_name]])
  draws <- seeded(SEED, replicate(N_RANDOM, sample(NON_TARGETS, n), simplify = FALSE))
  map_dfr(names(VARIANTS), function(v) {
    rhos <- vapply(draws, function(g) quick_rho(score_variant(g, v)), numeric(1))
    tibble(list = list_name, n_genes = n, variant = v,
           null_median = median(rhos), null_q95 = quantile(rhos, 0.95), null_q05 = quantile(rhos, 0.05))
  })
})
write.csv(null, file.path(OUT_DIR, "01_ladder_null.csv"), row.names = FALSE)

# ---- Summary printout: the induced-cell fluorescence column, side by side -----
cat("\n=== Induced cells, GFP fluorescence ===\n")
ladder |>
  filter(readout == "fluorescence", cell_set == "induced") |>
  select(list, variant, rho, ci_low, ci_high) |>
  left_join(select(null, list, variant, null_median, null_q95), by = c("list", "variant")) |>
  mutate(across(where(is.numeric), \(x) round(x, 3))) |>
  arrange(list, variant) |>
  print(n = Inf)

# ---- Figure -----------------------------------------------------------------
null_plot <- null |> mutate(readout = "fluorescence", cell_set = "induced")
p <- ladder |>
  mutate(variant = factor(variant, levels = names(VARIANTS)),
         panel = paste(readout, cell_set, sep = " | ")) |>
  ggplot(aes(x = variant, y = rho, colour = list, group = list)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey60") +
  geom_line(position = position_dodge(width = 0.3)) +
  geom_pointrange(aes(ymin = ci_low, ymax = ci_high), position = position_dodge(width = 0.3)) +
  geom_point(data = null_plot |> mutate(variant = factor(variant, levels = names(VARIANTS)),
                                        panel = "fluorescence | induced"),
             aes(y = null_q95, shape = "random genes, 95th pct"), position = position_dodge(width = 0.3)) +
  scale_shape_manual(values = c("random genes, 95th pct" = 4), name = NULL) +
  facet_wrap(~ panel) +
  labs(x = "Aggregation step (A = current method ... D = rank-based)", y = "Spearman rho with GFP",
       colour = "Target list",
       title = "Aggregation ladder: where does TarBase performance recover?",
       caption = paste(sprintf("%s = %s", names(VARIANTS), sub("^[A-Z0-9]+: ", "", VARIANTS)), collapse = "; ")) +
  theme_bw(base_size = 12)
ggsave(file.path(OUT_DIR, "01_ladder.pdf"), p, width = 12, height = 8)
cat("\nWrote 01_ladder_stratified.csv, 01_ladder_null.csv, 01_ladder.pdf to", OUT_DIR, "\n")
