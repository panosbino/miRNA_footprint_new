suppressMessages({ library(dplyr); library(tibble); library(tidyr); library(ggplot2) })

# ===========================================================================
# Accuracy of every method against GFP FLUORESCENCE (FACS index sort), on
# Dox-induced cells only.
#
# Three correlations are reported per method so the result can be interpreted:
#   A  vs GFP mRNA, all cells          -- the number reported before
#   B  vs GFP mRNA, induced cells      -- same cells as C; isolates the effect of
#                                         restricting to induced cells (A -> B)
#   C  vs GFP FLUORESCENCE, induced    -- the new, independent ground truth
#                                         (B -> C isolates the effect of the swap)
# The reference line is GFP mRNA itself against fluorescence: how well even the
# mRNA reporter tracks the protein signal.
#
# How cells are linked to fluorescence: FACS index-sort gives fluorescence per
# plate WELL; the sequencing barcode -> well map (P32156_<plate>.well_barcodes.txt)
# gives each cell's well. The join is validated empirically below (fluorescence
# must follow the Dox gradient across plate columns and correlate with GFP mRNA),
# because a wrong well geometry would scramble everything silently.
# ===========================================================================

BASE_DIR      <- Sys.getenv("HEK_BASE_DIR", "~/Desktop/Projects/miRNA_footprint_new")
RAW_DIR       <- file.path(BASE_DIR, "datasets/HEK_SS3/raw")
PROCESSED_DIR <- file.path(BASE_DIR, "datasets/HEK_SS3/processed")
COMP_DIR      <- file.path(BASE_DIR, "analysis/HEK_SS3/comparisons")
OUT_DIR       <- COMP_DIR
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

SEQ_PREFIX      <- "P32156"                 # prefix of the per-plate sequencing files (from the build script)
PLATES          <- c("101", "103", "105")
INDUCED_GROUPS  <- "Dox Induced"            # = 0.01, 0.1 and 1 Dox; excludes "Control", "0 Dox" and "Unknown"
GFP_UPPER_PCT   <- 0.99                     # the established GFP-mRNA trim, applied only to column A (as before)
USE_COMMON_CELLS <- TRUE                    # score every method on the SAME cells
MIN_MATCH_FRAC  <- 0.5                      # per plate: share of cells with a well that must find fluorescence
MIN_RHO_SANITY  <- 0.2                      # fluorescence vs GFP mRNA across all matched cells must exceed this
MIN_CELLS       <- 30

# --------------------------------------------------------------------------
# helpers
# --------------------------------------------------------------------------
# "A01", "A1", "a1" -> "A1"; anything that is not letters+digits -> NA
norm_well <- function(x) {
  x <- toupper(trimws(as.character(x)))
  ok <- grepl("^[A-Z]+[0-9]+$", x)
  out <- rep(NA_character_, length(x))
  out[ok] <- paste0(sub("^([A-Z]+)[0-9]+$", "\\1", x[ok]), as.integer(sub("^[A-Z]+([0-9]+)$", "\\1", x[ok])))
  out
}
# approximate 95% CI for a Spearman correlation (Bonett & Wright 2000)
spearman_ci <- function(r, n) { se <- sqrt((1 + r^2 / 2) / (n - 3)); tanh(atanh(r) + c(-1, 1) * 1.96 * se) }
sp <- function(x, y) suppressWarnings(cor(x, y, method = "spearman", use = "complete.obs"))

# --------------------------------------------------------------------------
# 1. cells, wells, fluorescence
# --------------------------------------------------------------------------
meta <- readRDS(file.path(PROCESSED_DIR, "cell_metadata.rds"))
need <- c("cell_id", "Dox_Group", "Dox_Concentration", "GFP_normalized")
if (length(setdiff(need, colnames(meta))) > 0) {
  stop("cell_metadata.rds lacks column(s): ", paste(setdiff(need, colnames(meta)), collapse = ", "),
       ". These are written by build_annotated_scran_dataset.R; this file may predate it.")
}
if (!"plate" %in% colnames(meta)) meta$plate <- sub("_.*", "", meta$cell_id)
meta <- meta %>% transmute(cell_id, plate = as.character(plate), Dox_Group = as.character(Dox_Group),
                           Dox_Concentration = as.character(Dox_Concentration), GFP_normalized = as.numeric(GFP_normalized))
cat(sprintf("cell_metadata: %d cells; Dox_Group: %s\n", nrow(meta),
            paste(names(table(meta$Dox_Group)), table(meta$Dox_Group), sep = "=", collapse = ", ")))

read_fluor <- function(plate) {
  f <- file.path(RAW_DIR, sprintf("P%s_GFP_fluorescence.tsv", plate))
  if (!file.exists(f)) { cat(sprintf("*** plate %s: %s not found -- this plate is EXCLUDED from the analysis ***\n", plate, basename(f))); return(NULL) }
  x <- read.delim(f, check.names = FALSE, stringsAsFactors = FALSE)
  if (!"GFP" %in% colnames(x)) stop(basename(f), " has no 'GFP' column; columns are: ", paste(colnames(x), collapse = ", "))
  if (!"Well_ID" %in% colnames(x)) {
    if (all(c("Row", "Column") %in% colnames(x))) x$Well_ID <- paste0(x$Row, x$Column)
    else stop(basename(f), " has neither 'Well_ID' nor 'Row'+'Column'; columns are: ", paste(colnames(x), collapse = ", "))
  }
  out <- tibble(plate = plate, well = norm_well(x$Well_ID), fluorescence = suppressWarnings(as.numeric(x$GFP)))
  n_bad <- sum(is.na(out$well)); if (n_bad > 0) cat(sprintf("plate %s: %d fluorescence rows with an unparseable well ID dropped\n", plate, n_bad))
  out <- out[!is.na(out$well), ]
  n_dup <- sum(duplicated(out$well)); if (n_dup > 0) cat(sprintf("plate %s: %d wells appear more than once in the fluorescence file; keeping the first\n", plate, n_dup))
  out[!duplicated(out$well), ]
}
read_wellmap <- function(plate) {
  f <- file.path(RAW_DIR, sprintf("%s_%s.well_barcodes.txt", SEQ_PREFIX, plate))
  if (!file.exists(f)) stop("Well map not found: ", f, ". Files containing 'well_barcodes' in ", RAW_DIR, ": ",
                            paste(list.files(RAW_DIR, pattern = "well_barcodes"), collapse = ", "))
  w <- read.delim(f, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  if (!all(c("bc_set", "WellID") %in% colnames(w))) stop(basename(f), " must have columns bc_set and WellID; found: ", paste(colnames(w), collapse = ", "))
  tibble(cell_id = paste0(plate, "_", w$bc_set), well = norm_well(w$WellID))
}

fl <- bind_rows(lapply(PLATES, read_fluor))
if (nrow(fl) == 0) stop("No fluorescence files found in ", RAW_DIR)
plates_used <- unique(fl$plate)
wm <- bind_rows(lapply(plates_used, read_wellmap))

cells <- meta %>% inner_join(wm, by = "cell_id") %>% left_join(fl, by = c("plate", "well"))

cat("\n--- Join diagnostics per plate ---\n")
diag <- cells %>% group_by(plate) %>%
  summarise(cells_with_well = n(), cells_with_fluorescence = sum(is.finite(fluorescence)), .groups = "drop") %>%
  mutate(match_frac = cells_with_fluorescence / cells_with_well)
print(as.data.frame(diag), row.names = FALSE)
bad <- diag$plate[diag$match_frac < MIN_MATCH_FRAC]
if (length(bad) > 0) {
  stop("Plate(s) ", paste(bad, collapse = ", "), ": fewer than ", 100 * MIN_MATCH_FRAC, "% of cells found a fluorescence value. ",
       "Well IDs may be formatted differently between the two files.\n  example well IDs in the well map: ",
       paste(head(unique(wm$well), 5), collapse = ", "), "\n  example well IDs in the fluorescence file: ", paste(head(unique(fl$well), 5), collapse = ", "))
}

# --------------------------------------------------------------------------
# 2. validate the join: the fluorescence must behave like the experiment
# --------------------------------------------------------------------------
matched <- cells %>% dplyr::filter(is.finite(fluorescence), is.finite(GFP_normalized))
rho_sanity <- sp(matched$fluorescence, matched$GFP_normalized)
dose_levels <- c("Control", "0 Dox", "0.01 Dox", "0.1 Dox", "1 Dox")
dose_med <- matched %>% dplyr::filter(Dox_Concentration %in% dose_levels) %>%
  group_by(Dox_Concentration) %>% summarise(n = n(), median_fluorescence = median(fluorescence), .groups = "drop") %>%
  arrange(match(Dox_Concentration, dose_levels))
cat(sprintf("\n--- Join validation ---\nSpearman(fluorescence, GFP mRNA), all matched cells (n = %d): %.3f\n", nrow(matched), rho_sanity))
cat("Median fluorescence by Dox concentration (should rise from 0 Dox to 1 Dox):\n"); print(as.data.frame(dose_med), row.names = FALSE)
if (!is.finite(rho_sanity) || rho_sanity < MIN_RHO_SANITY) {
  stop(sprintf("Fluorescence correlates only %.3f with GFP mRNA across matched cells (threshold %.1f). With the Dox gradient this should be strong, ",
               rho_sanity, MIN_RHO_SANITY), "so the well-to-cell join is almost certainly WRONG (e.g. rows/columns swapped between the FACS export and the ",
       "sequencing well map). Do not use any result from this script until that is resolved.")
}
dm <- dose_med$median_fluorescence[dose_med$Dox_Concentration %in% c("0 Dox", "0.01 Dox", "0.1 Dox", "1 Dox")]
if (length(dm) >= 3 && any(diff(dm) < 0)) {
  cat("*** WARNING: median fluorescence does not rise monotonically with Dox concentration -- check the join/plate layout before trusting results. ***\n")
} else cat("PASSED: fluorescence tracks both GFP mRNA and the Dox gradient.\n")

ind <- matched %>% dplyr::filter(Dox_Group %in% INDUCED_GROUPS)
cat(sprintf("\nInduced cells (%s) with a fluorescence value and GFP mRNA: %d\n", paste(INDUCED_GROUPS, collapse = ", "), nrow(ind)))
print(as.data.frame(ind %>% count(plate, Dox_Concentration) %>% pivot_wider(names_from = Dox_Concentration, values_from = n, values_fill = 0)), row.names = FALSE)
if (nrow(ind) < MIN_CELLS) stop("Too few induced cells with fluorescence (", nrow(ind), ").")

# --------------------------------------------------------------------------
# 3. per-cell scores for every method (saved result files)
# --------------------------------------------------------------------------
std_label <- "Our method (TargetScan top 200)"
specs <- list(
  list(file = "res_mireact_comparison_scran.rds",
       labels = c(our_targetscan = std_label, our_tarbase = "Our method (TarBase)", mireact = "miReact (TarBase mode)")),
  list(file = "res_bayesreact_comparison_scran.rds",
       labels = c(our = std_label, bayesreact = "bayesReact (motif)")),
  list(file = "res_mireact_default_comparison_scran.rds",
       labels = c(our = std_label, mireact_default = "miReact (motif mode)")),
  list(file = "res_mitea_comparison_scran.rds",
       labels = c(our = std_label, mitea_mirtarbase = "miTEA-HiRes (miRTarBase)", mitea_tarbase = "miTEA-HiRes (TarBase)",
                  mitea_targetscan_all = "miTEA-HiRes (TargetScan, all)", mitea_targetscan_top200 = "miTEA-HiRes (TargetScan, top 200)")))

entries <- list()
add_entry <- function(label, scores) {          # scores: named numeric vector (names = cell_id)
  scores <- scores[is.finite(scores)]
  if (label %in% names(entries)) {
    old <- entries[[label]]; common <- intersect(names(old), names(scores))
    if (length(common) > 0 && max(abs(old[common] - scores[common])) > 1e-6)
      cat(sprintf("*** WARNING: '%s' differs between result files; keeping the first. ***\n", label))
    return(invisible())
  }
  entries[[label]] <<- scores
}
cat("\n--- Loading saved per-cell scores ---\n")
for (sp_ in specs) {
  p <- file.path(COMP_DIR, sp_$file)
  if (!file.exists(p)) { cat(sprintf("SKIPPED (not found): %s\n", sp_$file)); next }
  r <- readRDS(p)
  for (el in names(sp_$labels)) {
    d <- r[[el]]
    if (!is.data.frame(d) || !"Cell_ID" %in% colnames(d)) { cat(sprintf("  %s: element '%s' is not a per-cell table; skipped\n", sp_$file, el)); next }
    num <- setdiff(colnames(d)[sapply(d, is.numeric)], "eGFP")
    if (length(num) == 0) { cat(sprintf("  %s/%s: no numeric score column; skipped\n", sp_$file, el)); next }
    if (length(num) > 1) cat(sprintf("  %s/%s: several numeric columns (%s); using '%s'\n", sp_$file, el, paste(num, collapse = ", "), num[1]))
    add_entry(sp_$labels[[el]], setNames(d[[num[1]]], d$Cell_ID))
  }
  other <- setdiff(names(r)[sapply(r, is.data.frame)], names(sp_$labels))
  if (length(other) > 0) cat(sprintf("  note: %s also holds per-cell tables not used here: %s\n", sp_$file, paste(other, collapse = ", ")))
}

# --------------------------------------------------------------------------
# 4. our method with the co-expression-selected targets (optional)
# --------------------------------------------------------------------------
cx_path <- file.path(COMP_DIR, "res_coexpression_validation.rds")
scran_path <- file.path(PROCESSED_DIR, "scran_normalized_linear.rds")
if (file.exists(cx_path) && file.exists(scran_path)) {
  cx <- readRDS(cx_path); counts <- as.matrix(readRDS(scran_path))
  # same score as calculate_activity(): -log2(sum of targets / mean over cells) + 1
  score_genes <- function(genes) {
    genes <- intersect(genes, rownames(counts)); s <- colSums(counts[genes, , drop = FALSE]); -log2(s / mean(s)) + 1
  }
  K <- length(cx$selected)
  add_entry(sprintf("Our method (co-expression selected, K=%d)", K), score_genes(cx$selected))
  if (is.data.frame(cx$features) && all(c("ensembl_gene_id", "context_score") %in% colnames(cx$features))) {
    ts_k <- cx$features$ensembl_gene_id[order(cx$features$context_score)][seq_len(K)]
    add_entry(sprintf("Our method (TargetScan top %d, same size)", K), score_genes(ts_k))
    add_entry(sprintf("Our method (all %d TargetScan candidates)", nrow(cx$features)), score_genes(cx$features$ensembl_gene_id))
  }
  cat(sprintf("Added co-expression variants (K = %d) computed from the scran matrix\n", K))
} else cat("Co-expression variants skipped (need res_coexpression_validation.rds and scran_normalized_linear.rds)\n")

if (length(entries) == 0) stop("No method scores found in ", COMP_DIR)

# --------------------------------------------------------------------------
# 5. accuracy table
# --------------------------------------------------------------------------
gfp_cut <- quantile(meta$GFP_normalized, GFP_UPPER_PCT, na.rm = TRUE)
common <- Reduce(intersect, lapply(entries, names))
analysis_cells <- if (USE_COMMON_CELLS) intersect(ind$cell_id, common) else ind$cell_id
cat(sprintf("\nCells scored by every method: %d. Induced cells with fluorescence: %d; analysed together: %d%s\n",
            length(common), nrow(ind), length(analysis_cells),
            ifelse(length(analysis_cells) < 0.95 * nrow(ind), "  (*** more than 5% of induced cells lost: methods exclude cells (GFP-mRNA trim / undefined scores) ***)", "")))
ind_a <- ind %>% dplyr::filter(cell_id %in% analysis_cells)
rho_ref <- sp(ind_a$fluorescence, ind_a$GFP_normalized)

res <- bind_rows(lapply(names(entries), function(lbl) {
  s <- entries[[lbl]]
  all_tbl <- tibble(cell_id = names(s), score = as.numeric(s)) %>% inner_join(meta[, c("cell_id", "GFP_normalized")], by = "cell_id") %>%
    dplyr::filter(GFP_normalized < gfp_cut)
  d <- tibble(cell_id = names(s), score = as.numeric(s)) %>%
    inner_join(if (USE_COMMON_CELLS) ind_a else ind, by = "cell_id")
  n <- nrow(d)
  rC <- sp(d$score, d$fluorescence); ci <- if (n > 3) spearman_ci(rC, n) else c(NA, NA)
  tibble(method = lbl, n_all_cells = nrow(all_tbl), rho_A_mRNA_all = sp(all_tbl$score, all_tbl$GFP_normalized),
         n_induced = n, rho_B_mRNA_induced = sp(d$score, d$GFP_normalized),
         rho_C_fluorescence = rC, C_ci_lo = ci[1], C_ci_hi = ci[2],
         p_C = if (n > 3) suppressWarnings(cor.test(d$score, d$fluorescence, method = "spearman", exact = FALSE)$p.value) else NA_real_)
})) %>% arrange(desc(rho_C_fluorescence))

cat(sprintf("\n=== Accuracy vs GFP fluorescence, induced cells (reference: GFP mRNA vs fluorescence on the same cells = %.3f) ===\n", rho_ref))
print(as.data.frame(res %>% mutate(across(c(starts_with("rho"), starts_with("C_ci")), ~ round(.x, 3)), p_C = signif(p_C, 2))), row.names = FALSE)
cat("\nA = vs GFP mRNA, all cells (as reported before) | B = vs GFP mRNA, induced cells | C = vs fluorescence, induced cells (95% CI approximate)\n",
    "A -> B isolates restricting to induced cells; B -> C isolates swapping mRNA for fluorescence. Differences of a few hundredths between\n",
    "methods are within the CIs and should not be read as a ranking.\n", sep = "")

# --------------------------------------------------------------------------
# 6. plots
# --------------------------------------------------------------------------
bar_df <- res %>% transmute(method, `GFP mRNA (same induced cells)` = rho_B_mRNA_induced, `GFP fluorescence (induced cells)` = rho_C_fluorescence) %>%
  pivot_longer(-method, names_to = "ground_truth", values_to = "rho") %>%
  mutate(method = factor(method, levels = rev(res$method)), ground_truth = factor(ground_truth, levels = c("GFP mRNA (same induced cells)", "GFP fluorescence (induced cells)")))
p_bar <- ggplot(bar_df, aes(method, rho, fill = ground_truth)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.75) +
  geom_text(aes(label = sprintf("%.2f", rho)), position = position_dodge(width = 0.8), hjust = -0.15, size = 3.3) +
  geom_hline(yintercept = rho_ref, linetype = "dashed") +
  coord_flip(clip = "off") + scale_fill_manual(values = c("grey70", "#2c7fb8")) +
  scale_y_continuous(expand = expansion(mult = c(0.02, 0.12))) +
  theme_bw(base_size = 13) + theme(legend.position = "bottom", legend.title = element_blank(), panel.grid.minor = element_blank()) +
  labs(x = NULL, y = "Spearman rho",
       title = "Method accuracy on Dox-induced cells",
       subtitle = sprintf("n = %d induced cells; ground truth: GFP mRNA (grey) vs GFP fluorescence (blue)", length(unique(ind_a$cell_id))),
       caption = sprintf("Dashed line: GFP mRNA vs GFP fluorescence on the same cells (rho = %.2f)", rho_ref))
ggsave(file.path(OUT_DIR, "fluorescence_accuracy_bars.pdf"), p_bar, width = 9.5, height = max(4.5, 0.5 * nrow(res) + 2.5))

chk <- matched %>% dplyr::filter(Dox_Concentration %in% dose_levels, fluorescence > 0) %>%
  mutate(Dox_Concentration = factor(Dox_Concentration, levels = dose_levels))
p_chk1 <- ggplot(chk, aes(Dox_Concentration, log10(fluorescence))) + geom_boxplot(fill = "grey85", outlier.alpha = 0.3) +
  theme_bw(base_size = 12) + labs(x = NULL, y = "log10 fluorescence (FITC-A)", title = "Join check: fluorescence by Dox concentration")
p_chk2 <- ggplot(chk, aes(log10(GFP_normalized + 1), log10(fluorescence), colour = Dox_Concentration)) + geom_point(alpha = 0.5, size = 1.2) +
  theme_bw(base_size = 12) + theme(legend.position = "bottom") +
  labs(x = "log10(GFP mRNA + 1)", y = "log10 fluorescence", title = sprintf("Join check: protein vs mRNA (Spearman %.2f)", rho_sanity))
pdf(file.path(OUT_DIR, "fluorescence_join_check.pdf"), width = 10, height = 5); print(p_chk1); print(p_chk2); invisible(dev.off())

write.csv(res, file.path(OUT_DIR, "fluorescence_accuracy_table.csv"), row.names = FALSE)
saveRDS(list(table = res, reference_rho = rho_ref, sanity_rho = rho_sanity, dose_medians = dose_med, join_diagnostics = diag,
             induced_cells = ind_a, induced_groups = INDUCED_GROUPS, plates_used = plates_used, use_common_cells = USE_COMMON_CELLS),
        file.path(OUT_DIR, "res_fluorescence_accuracy.rds"))
cat(sprintf("\nSaved table, plots and res_fluorescence_accuracy.rds to %s\n", OUT_DIR))

