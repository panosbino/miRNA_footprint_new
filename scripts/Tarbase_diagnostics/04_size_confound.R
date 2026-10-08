# ==============================================================================
# 04  Is the TarBase sum mostly measuring cell size / library size?
#
# (a) Covariates. Spearman rho of each score with log size factor, log total
#     UMIs and genes detected; Kruskal-Wallis p-value for plate. All and induced.
# (b) Generic-signal check. Score random non-target sets of TarBase's size with
#     the raw sum (A). If they correlate strongly with the TarBase score, the
#     TarBase sum mostly reflects a generic expression signal, not miR-124.
# (c) Residual. Regress each score on log size factor + log genes detected +
#     plate (all QC cells), then evaluate the residual against GFP.
#     CAUTION: if cell size is biologically linked to GFP level, this also
#     removes real signal. Read (c) together with (a), not on its own.
# (d) Fluorescence oddity. With TarBase, our score tracks fluorescence (0.27)
#     far better than GFP mRNA (0.06) in induced cells. FACS protein is not
#     normalized by cell size, so bigger cells carry more GFP. Checks whether
#     each readout tracks cell size, and the score-GFP correlation after removing
#     cell size from both (partial Spearman).
#
# Outputs: 04_covariates.csv, 04_random_sets.csv, 04_residualized.csv,
#          04_readout_covariates.csv, 04_partial.csv, 04_size_confound.pdf
# ==============================================================================

source(file.path(Sys.getenv("MIRNA_BASE_DIR", getwd()), "scripts/HEK_SS3/tarbase_diagnostics/tarbase_diag_common.R"))

SCORES <- list(
  `TarBase | A (raw sum)`           = score_variant(TARBASE, "A"),
  `TarBase | C (log ratio)`         = score_variant(TARBASE, "C"),
  `TarBase | D (rank)`              = score_variant(TARBASE, "D"),
  `TargetScan top200 | A (raw sum)` = score_variant(TS_TOP200, "A"),
  `TargetScan top200 | C (log ratio)` = score_variant(TS_TOP200, "C")
)

cells <- truth |> left_join(covariates, by = "Cell_ID") |>
  mutate(log_size_factor = log(size_factor), log_total_umis = log(total_umis),
         log_genes_detected = log(genes_detected))
COVS <- c("log_size_factor", "log_total_umis", "genes_detected")
subset_set <- function(d, cs) if (cs == "induced") d[d$induced, ] else d

# ---- (a) Score vs technical covariates ------------------------------------------
cov_res <- imap_dfr(SCORES, function(s, name) {
  map_dfr(CELL_SETS, function(cs) {
    d <- subset_set(cells, cs); d$score <- s[d$Cell_ID]
    bind_rows(
      map_dfr(COVS, \(cv) tibble(covariate = cv, rho = cor(d$score, d[[cv]], method = "spearman"), p_value = NA_real_)),
      tibble(covariate = "plate (Kruskal-Wallis p)", rho = NA_real_, p_value = kruskal.test(d$score, d$plate)$p.value)
    ) |> mutate(score = name, cell_set = cs)
  })
})
write.csv(cov_res, file.path(OUT_DIR, "04_covariates.csv"), row.names = FALSE)
cat("=== Score vs technical covariates (Spearman rho) ===\n")
cov_res |> filter(!is.na(rho)) |>
  select(score, cell_set, covariate, rho) |>
  pivot_wider(names_from = covariate, values_from = rho) |>
  mutate(across(where(is.double), \(x) round(x, 3))) |> print(n = Inf)

# ---- (b) Random sets of TarBase's size ------------------------------------------
cat(sprintf("\n=== %d random non-target sets of %d genes, raw sum (A) ===\n", N_RANDOM, length(TARBASE)))
draws <- seeded(SEED, replicate(N_RANDOM, sample(NON_TARGETS, min(length(TARBASE), length(NON_TARGETS))), simplify = FALSE))
rand <- map_dfr(seq_along(draws), function(i) {
  s <- score_variant(draws[[i]], "A")
  d <- cells; d$score <- s[d$Cell_ID]
  tibble(draw = i,
         rho_with_tarbase_score = cor(s[cells$Cell_ID], SCORES[[1]][cells$Cell_ID], method = "spearman"),
         rho_with_size_factor   = cor(d$score, d$log_size_factor, method = "spearman"),
         rho_facs_induced       = quick_rho(s, "gfp_facs", "induced"),
         rho_mrna_induced       = quick_rho(s, "gfp_mrna", "induced"))
})
write.csv(rand, file.path(OUT_DIR, "04_random_sets.csv"), row.names = FALSE)
rand |> summarise(across(-draw, list(median = median, min = min, max = max))) |>
  pivot_longer(everything()) |> mutate(value = round(value, 3)) |> print(n = Inf)

# ---- (c) Residualized scores ------------------------------------------------------
cat("\n=== Scores with cell size and plate regressed out ===\n")
resid_res <- imap_dfr(SCORES, function(s, name) {
  d <- cells; d$score <- s[d$Cell_ID]
  fit <- lm(score ~ log_size_factor + log_genes_detected + plate, data = d)
  r <- setNames(residuals(fit), d$Cell_ID)
  bind_rows(evaluate(s, paste(name, "| original"))  |> mutate(version = "original"),
            evaluate(r, paste(name, "| residual"))  |> mutate(version = "size + plate removed")) |>
    mutate(score = name, r2_covariates = summary(fit)$r.squared)
})
write.csv(resid_res, file.path(OUT_DIR, "04_residualized.csv"), row.names = FALSE)

# ---- (d) Fluorescence oddity --------------------------------------------------------
cat("\n=== Do the GFP readouts themselves track cell size? (Spearman rho) ===\n")
readout_cov <- map_dfr(CELL_SETS, function(cs) {
  d <- subset_set(cells, cs)
  expand_grid(readout = names(GFP_READOUTS), covariate = COVS) |>
    mutate(cell_set = cs,
           rho = map2_dbl(readout, covariate, \(ro, cv) cor(d[[GFP_READOUTS[[ro]]]], d[[cv]], method = "spearman")))
})
write.csv(readout_cov, file.path(OUT_DIR, "04_readout_covariates.csv"), row.names = FALSE)
readout_cov |> pivot_wider(names_from = covariate, values_from = rho) |>
  mutate(across(where(is.double), \(x) round(x, 3))) |> print()

# Partial Spearman: rank everything, remove covariate ranks from score and readout, correlate residuals.
partial_spearman <- function(x, y, Z) {
  rx <- rank(x); ry <- rank(y); RZ <- apply(as.matrix(Z), 2, rank)
  cor(residuals(lm(rx ~ RZ)), residuals(lm(ry ~ RZ)))
}
cat("\n=== Score vs GFP: plain vs partial Spearman (cell size removed), induced cells ===\n")
partial <- imap_dfr(SCORES, function(s, name) {
  d <- subset_set(cells, "induced"); d$score <- s[d$Cell_ID]
  map_dfr(names(GFP_READOUTS), function(ro) {
    y <- d[[GFP_READOUTS[[ro]]]]
    tibble(score = name, readout = ro,
           rho_plain   = cor(d$score, y, method = "spearman"),
           rho_partial = partial_spearman(d$score, y, d[, c("log_size_factor", "log_genes_detected")]))
  })
})
write.csv(partial, file.path(OUT_DIR, "04_partial.csv"), row.names = FALSE)
partial |> mutate(across(where(is.double), \(x) round(x, 3))) |> print(n = Inf)

# ---- Figure: original vs residualized --------------------------------------------
p <- resid_res |>
  mutate(panel = paste(readout, cell_set, sep = " | ")) |>
  ggplot(aes(y = score, x = rho, colour = version)) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey60") +
  geom_pointrange(aes(xmin = ci_low, xmax = ci_high), position = position_dodge(width = 0.5)) +
  facet_wrap(~ panel) +
  labs(x = "Spearman rho with GFP", y = NULL, colour = NULL,
       title = "Scores before and after removing cell size and plate",
       subtitle = "A large drop means the score was partly measuring cell size (or that size is linked to GFP)") +
  theme_bw(base_size = 11) + theme(legend.position = "bottom")
ggsave(file.path(OUT_DIR, "04_size_confound.pdf"), p, width = 12, height = 7)
cat("\nWrote 04_* files to", OUT_DIR, "\n")
