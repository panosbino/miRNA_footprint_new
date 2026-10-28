# ==============================================================================
# Step 7 (all datasets): benchmark metric for every method x target list, and
# its comparison with the negative control and the random lists.
#
#   Rscript 7_summarise.R <mESC|HEK|PSCSR>
#
# Ground truth and metric
#   mESC   KO vs Control, all cells           balanced error (%), lower = better
#   HEK    GFP fluorescence, all / induced     Spearman rho, higher = better
#   PSCSR  measured family expression, each    Spearman rho, higher = better
#          (cell line, family) pair within its line
#
# For each method, family, list type and evaluation:
#   benchmark          metric of the benchmarked list
#   negative_control   metric of the control miRNA's list
#   random_mean/_sd    over the N_RANDOM random lists
#   vs_control         benchmark minus control } oriented so that POSITIVE =
#   vs_random          benchmark minus random  } better than the control
# Methods without target lists: bayesReact and miReact motif mode (seed site;
# control = control site; random = random 7-mers), miDGD (predicted proportion of
# the family's miRNAs; control = control family; random = random model miRNAs).
#
# Outputs (<out_dir>/results/): results_all_lists.tsv, results.tsv, summary.tsv,
#   heatmaps.pdf, per_family_plots.pdf (one page per family and cell line / cell set)
# ==============================================================================
source(file.path(Sys.getenv("MIRNA_BASE_DIR", "/cfs/klemming/projects/supr/naiss2024-6-235/miRNA_footprint_new"),
                 "scripts/benchmark_v3/benchmark_config.R"))
source(file.path(BASE_DIR, "scripts/benchmark_v3/benchmark_lib.R"))
ds <- dataset_from_args()
res_dir <- file.path(ds$out_dir, "results"); dir.create(res_dir, showWarnings = FALSE)
set.seed(RANDOM_SEED)
tl <- readRDS(file.path(ds$out_dir, "target_lists.rds"))
NO_LIST <- "no list (motif / expression)"

read_if <- function(path, reader) {
  if (file.exists(path)) return(reader(path))
  cat("SKIPPED (not found):", path, "\n"); NULL
}

# ==============================================================================
# 1. Ground truth: one row per evaluation (which family, which cells, truth)
# ==============================================================================
if (ds$name == "mESC") {
  pheno <- read.csv(file.path(BASE_DIR, "datasets/mESCs/phenotype.csv"), fileEncoding = "UTF-8-BOM")
  u <- new.env(); sys.source(file.path(BASE_DIR, "scripts/mESCs/mESC_separation_utils.R"), envir = u)
  evaluations <- tibble(family = unique(tl$lists$family), eval_set = "all cells",
                        cells = list(pheno$cells), truth = list(pheno$Exp == "KO"))
  metric <- function(score, is_ko) 100 * u$sep_error(score, is_ko)$bal_error
  metric_name <- "balanced error (%)"; higher_is_better <- FALSE
} else if (ds$name == "HEK") {
  gfp <- read.csv(file.path(BASE_DIR, "datasets/HEK_SS3/processed/gfp_correlation_data.csv"))
  gfp <- gfp[is.finite(gfp$GFP_fluorescence), ]
  induced <- gfp$Dox_Group == "Dox Induced"
  evaluations <- tibble(family = unique(tl$lists$family), eval_set = c("all", "induced"),
                        cells = list(gfp$cell_id, gfp$cell_id[induced]),
                        truth = list(gfp$GFP_fluorescence, gfp$GFP_fluorescence[induced]))
  metric <- function(score, truth) cor(score, truth, method = "spearman")
  metric_name <- "Spearman rho vs GFP fluorescence"; higher_is_better <- TRUE
} else {
  data_dir <- file.path(ds$out_dir, "data")
  cell_table <- read_tsv(file.path(data_dir, "cells.tsv"), show_col_types = FALSE)
  selection <- read_tsv(file.path(data_dir, "selection.tsv"), show_col_types = FALSE)
  family_expr <- readRDS(file.path(data_dir, "family_expression.rds"))
  line_cells <- split(cell_table$cell, cell_table$cell_line)
  evaluations <- tibble(family = selection$family, eval_set = selection$cell_line,
                        cells = line_cells[selection$cell_line],
                        truth = map2(selection$family, selection$cell_line, \(f, l) family_expr[f, line_cells[[l]]]))
  metric <- function(score, truth) cor(score, truth, method = "spearman")
  metric_name <- "Spearman rho vs measured miRNA (within line)"; higher_is_better <- TRUE
}
cat(sprintf("%s: %d evaluations; metric: %s\n", ds$name, nrow(evaluations), metric_name))

# ==============================================================================
# 2. Scores: one matrix per method (rows = list rows, columns = cells) plus a
#    table describing each row (family, list_type, kind, draw)
# ==============================================================================
score_dir <- file.path(ds$out_dir, "scores")
list_info <- tl$lists |> dplyr::select(list_id, family, list_type, kind, draw, n_genes)
methods <- list()                                   # name -> list(matrix, info)

for (m in c(paste0("Ours-", c("A", "B", "C", "C2", "D")), "miReact")) {
  S <- read_if(file.path(score_dir, paste0(m, ".rds")), readRDS)
  if (!is.null(S)) methods[[m]] <- list(scores = S, info = list_info)
}
mt <- read_if(file.path(score_dir, "miTEA-HiRes.tsv"), \(f) data.table::fread(f, data.table = FALSE))
if (!is.null(mt)) {
  S <- as.matrix(mt[, -1]); rownames(S) <- mt$list_id
  methods[["miTEA-HiRes"]] <- list(scores = S, info = filter(list_info, list_id %in% rownames(S)))
}

# Motif-based methods: a row's score = sum over its motifs (mESC benchmarks three sites together)
motif_rows <- tl$motif_lists |> mutate(list_id = sprintf("M%04d", row_number()), list_type = NO_LIST, n_genes = NA_integer_)
sum_motifs <- function(get_activity) {
  all_motifs <- unique(unlist(motif_rows$motifs))
  act <- map(set_names(all_motifs), get_activity)
  missing <- names(act)[map_lgl(act, is.null)]
  if (length(missing) == length(all_motifs)) return(NULL)
  if (length(missing)) cat("  motifs without results (rows using them are NA):", paste(missing, collapse = ", "), "\n")
  cells_all <- names(act[[setdiff(names(act), missing)[1]]])
  S <- t(vapply(motif_rows$motifs, function(ms) {
    if (any(ms %in% missing)) return(rep(NA_real_, length(cells_all)))
    Reduce(`+`, map(ms, \(x) act[[x]][cells_all]))
  }, numeric(length(cells_all))))
  dimnames(S) <- list(motif_rows$list_id, cells_all); S
}
cat("bayesReact:\n")
S <- sum_motifs(function(motif) {
  r <- read_if(file.path(ds$out_dir, "bayesreact", paste0(motif, ".rds")), readRDS)
  if (is.null(r)) return(NULL)
  if (max(r$diagnostics$Rhat, na.rm = TRUE) > 1.01) cat(sprintf("  *** %s: max Rhat %.3f (convergence doubtful)\n", motif, max(r$diagnostics$Rhat, na.rm = TRUE)))
  r$activity
})
if (!is.null(S)) methods[["bayesReact"]] <- list(scores = S, info = dplyr::select(motif_rows, list_id, family, list_type, kind, draw, n_genes))
cat("miReact (motif):\n")
job <- read_if(file.path(ds$out_dir, "mireact_motif_job.rds"), readRDS)
motif_matrix <- if (!is.null(job)) read_if(file.path(job$wd, job$out_file), readRDS)
if (!is.null(motif_matrix)) {
  S <- sum_motifs(\(motif) if (motif %in% rownames(motif_matrix)) motif_matrix[motif, ] else NULL)
  methods[["miReact (motif)"]] <- list(scores = S, info = dplyr::select(motif_rows, list_id, family, list_type, kind, draw, n_genes))
}

# miDGD (human only): predicted proportion of the family's miRNAs
midgd_file <- file.path(ds$out_dir, "midgd", "4a_mirna_prop.tsv")
if (DATASETS[[ds$name]]$species == "hs") {
  mg <- read_if(midgd_file, \(f) data.table::fread(f, data.table = FALSE))
  if (!is.null(mg)) {
    P <- as.matrix(mg[, -1]); rownames(P) <- mg[[1]]
    fam <- tl$families
    ctrl_members <- intersect(unlist(fam$members[fam$role == "control"]), colnames(P))
    rows <- list(); info <- list()
    for (f in unique(tl$lists$family)) {
      own <- intersect(unlist(fam$members[fam$role == "benchmark" & fam$family == f]), colnames(P))
      if (!length(own)) { cat("  miDGD: no model miRNA for", f, "\n"); next }
      random_mirnas <- sample(setdiff(colnames(P), c(own, ctrl_members)), N_RANDOM)
      cols <- c(list(own, ctrl_members), as.list(random_mirnas))
      for (k in seq_along(cols)) rows[[length(rows) + 1]] <- if (length(cols[[k]])) rowSums(P[, cols[[k]], drop = FALSE]) else rep(NA_real_, nrow(P))
      info[[f]] <- tibble(family = f, kind = c("benchmark", "negative_control", rep("random", N_RANDOM)),
                          draw = c(NA, NA, seq_len(N_RANDOM)))
    }
    if (length(rows)) {
      info <- bind_rows(info) |> mutate(list_id = sprintf("D%04d", row_number()), list_type = NO_LIST, n_genes = NA_integer_)
      S <- do.call(rbind, rows); dimnames(S) <- list(info$list_id, rownames(P))
      methods[["miDGD (4a)"]] <- list(scores = S, info = info)
      if (!length(ctrl_members)) cat("  miDGD: the control family is not in the model -- no negative control for miDGD\n")
    }
  }
}
cat(sprintf("Methods with scores: %s\n", paste(names(methods), collapse = ", ")))

# ==============================================================================
# 3. Metric for every method x list row x evaluation of that row's family
# ==============================================================================
safe_metric <- function(score, truth) {
  ok <- is.finite(score)
  if (sum(ok) < 10 || length(unique(score[ok])) < 2) return(NA_real_)
  metric(score[ok], truth[ok])
}
results_all <- imap_dfr(methods, function(m, method) {
  map_dfr(seq_len(nrow(evaluations)), function(e) {
    ev <- evaluations[e, ]
    rows <- m$info |> filter(family == ev$family)
    cc <- intersect(ev$cells[[1]], colnames(m$scores)); truth <- ev$truth[[1]][match(cc, ev$cells[[1]])]
    rows |> mutate(method = method, eval_set = ev$eval_set,
                   value = map_dbl(list_id, \(id) safe_metric(m$scores[id, cc], truth)))
  })
})
write_tsv(results_all, file.path(res_dir, "results_all_lists.tsv"))

# ==============================================================================
# 4. Benchmark vs negative control vs random
# ==============================================================================
better_by <- function(x, reference) if (higher_is_better) x - reference else reference - x
results <- results_all |>
  group_by(method, family, list_type, eval_set) |>
  summarise(benchmark        = value[kind == "benchmark"][1],
            negative_control = value[kind == "negative_control"][1],
            random_mean      = if (any(kind == "random")) mean(value[kind == "random"], na.rm = TRUE) else NA_real_,
            random_sd        = if (any(kind == "random")) sd(value[kind == "random"], na.rm = TRUE) else NA_real_,
            n_random         = sum(kind == "random" & is.finite(value)),
            n_genes_control  = n_genes[kind == "negative_control"][1],   # before n_genes, which summarise() overwrites
            n_genes          = n_genes[kind == "benchmark"][1],
            .groups = "drop") |>
  mutate(vs_control = better_by(benchmark, negative_control),
         vs_random  = better_by(benchmark, random_mean),
         rankable   = list_type == NO_LIST | (n_genes >= MIN_LIST_SIZE & !is.na(n_genes)))
write_tsv(results, file.path(res_dir, "results.tsv"))

# Summary: median over families and evaluations. Kept separate per cell set for HEK,
# and per benchmark for mESC (combined set and the three single families).
group_cols <- switch(ds$name, HEK = c("method", "list_type", "eval_set"),
                     mESC = c("method", "list_type", "family"), c("method", "list_type"))
summary_tbl <- results |> filter(rankable) |>
  group_by(across(all_of(group_cols))) |>
  summarise(n = n(),
            # shares first: summarise() replaces vs_control / vs_random with their medians below
            share_beating_control = mean(vs_control > 0, na.rm = TRUE),
            share_beating_random  = mean(vs_random > 0, na.rm = TRUE),
            benchmark = median(benchmark, na.rm = TRUE),
            negative_control = median(negative_control, na.rm = TRUE), random_mean = median(random_mean, na.rm = TRUE),
            vs_control = median(vs_control, na.rm = TRUE), vs_random = median(vs_random, na.rm = TRUE), .groups = "drop") |>
  arrange(desc(vs_random))
write_tsv(summary_tbl, file.path(res_dir, "summary.tsv"))
cat(sprintf("\n=== %s: %s. vs_control / vs_random > 0 = better than the control ===\n", ds$name, metric_name))
print(summary_tbl |> mutate(across(where(is.numeric), \(x) round(x, 3))), n = 60, width = Inf)

# ==============================================================================
# 5. Heatmaps: raw metric, difference to negative control, difference to random
# ==============================================================================
method_levels <- c(paste0("Ours-", c("A", "B", "C", "C2", "D")), "miReact", "miTEA-HiRes", "miReact (motif)", "bayesReact", "miDGD (4a)")
type_levels <- c(unique(tl$lists$list_type), NO_LIST)
heatmap <- function(col, title, diverging) {
  d <- summary_tbl |> mutate(method = factor(method, rev(method_levels)), list_type = factor(list_type, type_levels))
  if ("family" %in% names(d)) d <- mutate(d, family = factor(family, unique(tl$lists$family)))   # configured order: combined first
  p <- ggplot(d, aes(list_type, method, fill = .data[[col]])) +
    geom_tile(colour = "white") +
    geom_text(aes(label = ifelse(is.na(.data[[col]]), "", sprintf(if (higher_is_better) "%.2f" else "%.1f", .data[[col]]))), size = 2.4) +
    theme_bw(base_size = 10) + theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(x = NULL, y = NULL, fill = NULL, title = title,
         subtitle = sprintf("%s; %s; lists < %d genes excluded", metric_name,
                            switch(ds$name, PSCSR = "median over (cell line, family) pairs", mESC = "one panel per benchmark",
                                   "one panel per cell set"), MIN_LIST_SIZE))
  p <- p + if (diverging) scale_fill_gradient2(low = "#b2182b", mid = "white", high = "#2166ac", midpoint = 0, na.value = "grey90")
           else scale_fill_viridis_c(direction = if (higher_is_better) 1 else -1, na.value = "grey90")
  if (ds$name == "HEK") return(p + facet_wrap(~ eval_set, ncol = 1))
  if (ds$name == "mESC") return(p + facet_wrap(~ family, ncol = 1))
  p
}
pdf(file.path(res_dir, "heatmaps.pdf"), width = 18, height = switch(ds$name, HEK = 10, mESC = 18, 6))
print(heatmap("benchmark", "Benchmark metric (raw)", diverging = FALSE))
print(heatmap("vs_control", "Difference to the negative control (positive = better than control)", diverging = TRUE))
print(heatmap("vs_random", "Difference to random lists (positive = better than random)", diverging = TRUE))
invisible(dev.off())
# ==============================================================================
# 6. One page per family and evaluation (PSCSR: cell line; HEK: cell set):
#    every method and list type, with both controls next to the benchmarked value
# ==============================================================================
chance <- if (higher_is_better) 0 else 50
plot_family <- function(fam, ev_set) {
  d <- results |> filter(family == fam, eval_set == ev_set)
  if (!nrow(d)) return(NULL)
  d <- d |> mutate(method = factor(method, method_levels), list_type = factor(list_type, type_levels),
                   size_class = ifelse(rankable, "benchmarked list", sprintf("benchmarked list (< %d genes)", MIN_LIST_SIZE)))
  ggplot(d, aes(x = list_type)) +
    geom_hline(yintercept = chance, colour = "grey60", linetype = 2) +
    geom_linerange(aes(ymin = random_mean - random_sd, ymax = random_mean + random_sd),
                   colour = "grey60", linewidth = 2.2, alpha = 0.6) +
    geom_point(aes(y = random_mean), colour = "grey30", shape = 95, size = 4) +
    geom_point(aes(y = negative_control), colour = "#d6604d", shape = 4, size = 1.8, stroke = 0.9) +
    geom_point(aes(y = benchmark, shape = size_class), colour = "#2166ac", size = 2) +
    scale_shape_manual(values = setNames(c(16, 1), c("benchmarked list", sprintf("benchmarked list (< %d genes)", MIN_LIST_SIZE))),
                       name = NULL) +
    facet_wrap(~ method, ncol = 2) +
    theme_bw(base_size = 9) +
    theme(axis.text.x = element_text(angle = 60, hjust = 1), legend.position = "top") +
    labs(x = NULL, y = metric_name, title = sprintf("%s  |  %s", fam, ev_set),
         subtitle = "Blue = benchmarked list; red x = negative control; grey bar = random lists (mean +/- SD); dashed = chance")
}
pages <- distinct(evaluations, family, eval_set)                    # PSCSR: cell line, then rank, as chosen in step 0
pdf(file.path(res_dir, "per_family_plots.pdf"), width = 12, height = 11)
for (k in seq_len(nrow(pages))) {
  p <- plot_family(pages$family[k], pages$eval_set[k])
  if (!is.null(p)) suppressWarnings(print(p))                       # missing controls (NA) are simply not drawn
}
invisible(dev.off())

cat(sprintf("\nWrote results_all_lists.tsv, results.tsv, summary.tsv, heatmaps.pdf, per_family_plots.pdf to %s\n", res_dir))
