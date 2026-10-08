# ==============================================================================
# 02  Weighting: do a few highly expressed genes dominate the raw sum?
#
# (a) Concentration. For each cell: share of the target sum coming from that
#     cell's top 10 and top 50 contributing genes. Plus the list-level curve:
#     genes sorted by mean expression vs cumulative share of the summed mean.
# (b) Drop the top. Remove the 1, 5, 10, 25% most highly expressed genes of each
#     list and recompute the raw sum (A) and the log-ratio score (C).
#
# Reading the result: if TarBase's sum is driven by a handful of genes AND
# removing them raises rho for A (but barely moves C), weighting is the cause.
# If removing them lowers rho, those genes carry the signal instead.
#
# Outputs:
#   02_concentration_per_cell.csv, 02_concentration.pdf
#   02_drop_top_genes.csv,         02_drop_top_genes.pdf
# ==============================================================================

source(file.path(Sys.getenv("MIRNA_BASE_DIR", getwd()), "scripts/HEK_SS3/tarbase_diagnostics/tarbase_diag_common.R"))

# ---- (a) Concentration --------------------------------------------------------
top_share <- function(genes, k) {
  M <- X[genes, , drop = FALSE]
  apply(M, 2, function(x) { tot <- sum(x); if (tot == 0) NA_real_ else sum(sort(x, decreasing = TRUE)[seq_len(min(k, length(x)))]) / tot })
}
conc <- map_dfr(names(LISTS), function(list_name) {
  g <- LISTS[[list_name]]
  tibble(list = list_name, n_genes = length(g), Cell_ID = colnames(X),
         share_top10 = top_share(g, 10), share_top50 = top_share(g, 50))
})
write.csv(conc, file.path(OUT_DIR, "02_concentration_per_cell.csv"), row.names = FALSE)

cat("=== Share of each cell's target sum from its top-k genes (median [IQR] across cells) ===\n")
conc |>
  group_by(list, n_genes) |>
  summarise(across(c(share_top10, share_top50),
                   \(x) sprintf("%.2f [%.2f-%.2f]", median(x, na.rm = TRUE),
                                quantile(x, .25, na.rm = TRUE), quantile(x, .75, na.rm = TRUE))),
            .groups = "drop") |>
  print()

cum_curve <- map_dfr(names(LISTS), function(list_name) {
  m <- sort(gene_mean[LISTS[[list_name]]], decreasing = TRUE)
  tibble(list = list_name, gene_fraction = seq_along(m) / length(m), cumulative_share = cumsum(m) / sum(m))
})
p_conc <- ggplot(cum_curve, aes(gene_fraction, cumulative_share, colour = list)) +
  geom_abline(linetype = "dashed", colour = "grey60") +
  geom_line(linewidth = 1) +
  labs(x = "Fraction of list genes (most highly expressed first)",
       y = "Cumulative share of the summed mean expression",
       title = "How concentrated is each list's raw sum?",
       subtitle = "Dashed line = every gene contributes equally") +
  theme_bw(base_size = 12)
ggsave(file.path(OUT_DIR, "02_concentration.pdf"), p_conc, width = 8, height = 6)

# ---- (b) Drop the most highly expressed genes ---------------------------------
DROP_FRACTIONS <- c(0, 0.01, 0.05, 0.10, 0.25)
cat("\n=== Dropping the most highly expressed genes (scores A and C) ===\n")
drop_res <- map_dfr(names(LISTS), function(list_name) {
  g_sorted <- names(sort(gene_mean[LISTS[[list_name]]], decreasing = TRUE))
  map_dfr(DROP_FRACTIONS, function(f) {
    g <- g_sorted[(floor(f * length(g_sorted)) + 1):length(g_sorted)]
    map_dfr(c("A", "C"), function(v)
      evaluate(score_variant(g, v), sprintf("%s | %s | drop top %.0f%%", list_name, v, 100 * f)) |>
        mutate(list = list_name, variant = v, drop_fraction = f, n_genes = length(g)))
  })
})
write.csv(drop_res, file.path(OUT_DIR, "02_drop_top_genes.csv"), row.names = FALSE)

p_drop <- drop_res |>
  mutate(panel = paste(readout, cell_set, sep = " | "), variant = VARIANTS[variant]) |>
  ggplot(aes(100 * drop_fraction, rho, colour = list, linetype = variant)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey60") +
  geom_line() +
  geom_pointrange(aes(ymin = ci_low, ymax = ci_high), size = 0.3) +
  facet_wrap(~ panel) +
  labs(x = "% of most highly expressed list genes removed", y = "Spearman rho with GFP",
       colour = "Target list", linetype = "Score",
       title = "Does removing highly expressed genes help the raw sum?") +
  theme_bw(base_size = 12)
ggsave(file.path(OUT_DIR, "02_drop_top_genes.pdf"), p_drop, width = 12, height = 8)
cat("\nWrote 02_* files to", OUT_DIR, "\n")
