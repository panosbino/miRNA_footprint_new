# ==============================================================================
# 03  List quality: how many TarBase genes respond to miR-124 in these cells?
#
# (a) Per-gene response. Each gene's own Spearman rho with GFP fluorescence in
#     induced cells, for TarBase, TargetScan top 200 and random non-targets.
#     Responsive targets should go DOWN as miR-124 goes up (negative rho).
# (b) Overlap. TarBase genes that also have a TargetScan site vs TarBase-only.
# (c) Evidence. TarBase genes supported by CLIP-type, microarray, RNA-seq or
#     proteomics experiments, and by HEK293 vs other cell lines. Each subset is
#     compared with random TarBase subsets of the same size.
# (d) List size. Random TarBase subsets of increasing size vs the ranked
#     TargetScan list at the same sizes, scored with A (sum), C (log ratio),
#     D (rank). If A degrades with size while C/D don't, sum aggregation scales
#     badly with long lists.
#
# Outputs: 03_per_gene_rho.csv/.pdf, 03_subsets.csv/.pdf, 03_size_curve.csv/.pdf
# ==============================================================================

source(file.path(Sys.getenv("MIRNA_BASE_DIR", getwd()), "scripts/HEK_SS3/tarbase_diagnostics/tarbase_diag_common.R"))

induced_cells <- truth$Cell_ID[truth$induced]
facs_induced  <- truth$gfp_facs[truth$induced]

# ---- (a) Per-gene rho -----------------------------------------------------------
per_gene_rho <- function(genes) {
  M <- t(X[genes, induced_cells, drop = FALSE])
  r <- suppressWarnings(cor(M, facs_induced, method = "spearman"))[, 1]
  tibble(gene = genes, rho = r, detection = colMeans(M > 0), mean_expr = gene_mean[genes])
}
random_genes <- seeded(SEED, sample(NON_TARGETS, min(2000, length(NON_TARGETS))))
pg <- bind_rows(
  per_gene_rho(TARBASE)      |> mutate(set = "TarBase"),
  per_gene_rho(TS_TOP200)    |> mutate(set = "TargetScan top200"),
  per_gene_rho(random_genes) |> mutate(set = "Random non-targets")
)
write.csv(pg, file.path(OUT_DIR, "03_per_gene_rho.csv"), row.names = FALSE)

null_q05 <- quantile(pg$rho[pg$set == "Random non-targets"], 0.05, na.rm = TRUE)
cat("=== Per-gene rho with fluorescence, induced cells ===\n")
cat(sprintf("Threshold for 'responsive': below the 5th percentile of random genes (rho < %.3f)\n", null_q05))
pg |>
  group_by(set) |>
  summarise(n_genes = n(), n_tested = sum(!is.na(rho)),
            median_rho = median(rho, na.rm = TRUE),
            frac_responsive = mean(rho < null_q05, na.rm = TRUE),
            median_detection = median(detection), .groups = "drop") |>
  mutate(across(where(is.double), \(x) round(x, 3))) |>
  print()

p_pg <- ggplot(filter(pg, !is.na(rho)), aes(rho, colour = set)) +
  geom_density(linewidth = 1) +
  geom_vline(xintercept = c(0, null_q05), linetype = c("dashed", "dotted"), colour = "grey50") +
  labs(x = "Per-gene Spearman rho with GFP fluorescence (induced cells)", y = "Density", colour = NULL,
       title = "Do list genes go down when miR-124 goes up?",
       subtitle = "Dotted line = 5th percentile of random non-target genes") +
  theme_bw(base_size = 12)
ggsave(file.path(OUT_DIR, "03_per_gene_rho.pdf"), p_pg, width = 9, height = 6)

# Mean expression vs per-gene rho: are the genes that dominate the sum responsive?
p_pg2 <- ggplot(filter(pg, set == "TarBase", !is.na(rho)), aes(log10(mean_expr), rho)) +
  geom_point(alpha = 0.3, size = 0.8) + geom_smooth(method = "loess", se = TRUE, colour = "firebrick") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(x = "log10 mean expression", y = "Per-gene rho with fluorescence (induced)",
       title = "TarBase genes: response vs expression level",
       subtitle = "Highly expressed genes dominate the raw sum; are they responsive?") +
  theme_bw(base_size = 12)
ggsave(file.path(OUT_DIR, "03_per_gene_rho_vs_expression.pdf"), p_pg2, width = 8, height = 6)

# ---- (b)+(c) Subsets --------------------------------------------------------------
ev <- tarbase_entries[tarbase_entries$gid %in% TARBASE, ]
method_class <- case_when(
  grepl("CLIP|CLASH|AGO|IMPACT", ev$method, ignore.case = TRUE) ~ "CLIP / AGO-IP",
  grepl("Microarray", ev$method, ignore.case = TRUE)            ~ "Microarray",
  grepl("RNA-Seq", ev$method, ignore.case = TRUE)               ~ "RNA-seq",
  grepl("SILAC|proteom|Western", ev$method, ignore.case = TRUE) ~ "Proteomics / protein",
  TRUE                                                          ~ "Other")
is_hek <- grepl("HEK|293", ev$cell_line, ignore.case = TRUE)

overlap_label <- if (TS_FULL_LIST) "any TargetScan site" else "TargetScan top 1000"
SUBSETS <- c(
  setNames(list(intersect(TARBASE, TS_ANY), setdiff(TARBASE, TS_ANY)),
           c(paste("TarBase with", overlap_label), paste("TarBase without", overlap_label))),
  split(ev$gid, paste("Evidence:", method_class)) |> lapply(unique),
  list(`Cell line: any HEK293 entry` = unique(ev$gid[is_hek]),
       `Cell line: no HEK293 entry`  = setdiff(TARBASE, ev$gid[is_hek]))
)
SUBSETS <- SUBSETS[lengths(SUBSETS) >= 20]
cat("\n=== TarBase subsets (genes are in several evidence subsets when supported by several experiments) ===\n")
print(tibble(subset = names(SUBSETS), n_genes = lengths(SUBSETS)))

subset_res <- imap_dfr(SUBSETS, function(g, name) {
  map_dfr(c("A", "C"), function(v) {
    obs <- evaluate(score_variant(g, v), sprintf("%s | %s", name, v))
    # Reference: random TarBase subsets of the same size (fluorescence, induced only)
    draws <- seeded(SEED, replicate(N_RANDOM, sample(TARBASE, length(g)), simplify = FALSE))
    ref <- vapply(draws, function(d) quick_rho(score_variant(d, v)), numeric(1))
    obs |> mutate(subset = name, variant = v, n_genes = length(g),
                  same_size_tarbase_median = median(ref),
                  same_size_tarbase_q05 = quantile(ref, .05), same_size_tarbase_q95 = quantile(ref, .95))
  })
})
write.csv(subset_res, file.path(OUT_DIR, "03_subsets.csv"), row.names = FALSE)

p_sub <- subset_res |>
  filter(readout == "fluorescence", cell_set == "induced") |>
  mutate(label = sprintf("%s (n=%d)", subset, n_genes), variant = VARIANTS[variant]) |>
  ggplot(aes(y = label, colour = variant)) +
  geom_linerange(aes(xmin = same_size_tarbase_q05, xmax = same_size_tarbase_q95),
                 linewidth = 4, alpha = 0.25, position = position_dodge(width = 0.6)) +
  geom_pointrange(aes(x = rho, xmin = ci_low, xmax = ci_high), position = position_dodge(width = 0.6)) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey60") +
  labs(x = "Spearman rho with fluorescence (induced cells)", y = NULL, colour = "Score",
       title = "TarBase subsets vs random TarBase subsets of the same size",
       subtitle = "Point + line = subset with 95% CI; shaded bar = 5th-95th pct of same-size random TarBase subsets") +
  theme_bw(base_size = 11)
ggsave(file.path(OUT_DIR, "03_subsets.pdf"), p_sub, width = 11, height = 0.6 * length(SUBSETS) + 2.5)

# ---- (d) Size curve ------------------------------------------------------------------
SIZES <- c(50, 100, 200, 500, 1000, 2000)
cat(sprintf("\n=== Size curve: %d random TarBase subsets per size; ranked TargetScan at the same sizes ===\n", N_RANDOM))
size_curve <- map_dfr(c("A", "C", "D"), function(v) {
  tb <- map_dfr(c(SIZES[SIZES < length(TARBASE)], length(TARBASE)), function(n) {
    draws <- if (n == length(TARBASE)) list(TARBASE) else
      seeded(SEED + n, replicate(N_RANDOM, sample(TARBASE, n), simplify = FALSE))
    r <- vapply(draws, function(d) quick_rho(score_variant(d, v)), numeric(1))
    tibble(list = "TarBase (random subsets)", size = n, rho_median = median(r),
           rho_q25 = quantile(r, .25), rho_q75 = quantile(r, .75))
  })
  ts_sizes <- c(SIZES[SIZES < length(TS_ANY)], length(TS_ANY))
  ts <- map_dfr(ts_sizes, function(n) {
    r <- quick_rho(score_variant(intersect(TS_RANKED, rownames(X))[1:n], v))
    tibble(list = "TargetScan (top-N, ranked)", size = n, rho_median = r, rho_q25 = r, rho_q75 = r)
  })
  bind_rows(tb, ts) |> mutate(variant = v)
})
write.csv(size_curve, file.path(OUT_DIR, "03_size_curve.csv"), row.names = FALSE)
size_curve |> mutate(across(where(is.double), \(x) round(x, 3))) |> print(n = Inf)

p_size <- size_curve |>
  mutate(variant = VARIANTS[variant]) |>
  ggplot(aes(size, rho_median, colour = variant, linetype = list)) +
  geom_ribbon(aes(ymin = rho_q25, ymax = rho_q75, fill = variant), alpha = 0.15, colour = NA) +
  geom_line() + geom_point() +
  scale_x_log10() +
  labs(x = "Number of target genes (log scale)", y = "Spearman rho with fluorescence (induced cells)",
       colour = "Score", fill = "Score", linetype = NULL,
       title = "Does the raw sum degrade with list size when log ratio / rank don't?",
       subtitle = "TarBase: median and IQR over random subsets") +
  theme_bw(base_size = 12)
ggsave(file.path(OUT_DIR, "03_size_curve.pdf"), p_size, width = 10, height = 6)
cat("\nWrote 03_* files to", OUT_DIR, "\n")
