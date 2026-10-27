# ==============================================================================
# PSCSR-seq V2: per-cell QC by cell line, for the mRNA and the miRNA library.
#
# Per cell (raw counts from datasets/SC_PSCSR_v2):
#   total_counts     sum of counts (UMIs; mRNA has multi-mapping weights)
#   detected         features with count > 0
#   zero             features with count 0 (features in the table minus detected)
# Per cell line: number of cells, median / mean / IQR of each, plus features
#   never detected in any cell of the line and features detected in every cell.
#
# Outputs (analysis/PSCSR/qc/): qc_per_cell.tsv, qc_per_cell_line.tsv, qc_per_cell_line.pdf
# ==============================================================================
source(file.path(Sys.getenv("MIRNA_BASE_DIR", "/cfs/klemming/projects/supr/naiss2024-6-235/miRNA_footprint_new"),
                 "scripts/benchmark_v2/benchmark_config.R"))
source(file.path(BASE_DIR, "scripts/benchmark_v2/benchmark_lib.R"))
out_dir <- file.path(BASE_DIR, "analysis/PSCSR/qc"); dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
inp <- function(f) file.path(PSCSR_INPUT$dir, PSCSR_INPUT[[f]])

anno  <- data.table::fread(inp("cell_anno"), data.table = FALSE)
cells <- tibble(cell = as.character(anno[[PSCSR_INPUT$cell_col]]), cell_line = as.character(anno[[PSCSR_INPUT$line_col]]))
libraries <- list(mRNA = read_count_matrix(inp("mrna_counts")), miRNA = read_count_matrix(inp("mirna_counts")))

# ---- Per cell ------------------------------------------------------------------------------------
per_cell <- imap_dfr(libraries, function(m, lib) {
  m <- m[, intersect(cells$cell, colnames(m)), drop = FALSE]
  tibble(library = lib, cell = colnames(m), features_in_table = nrow(m),
         total_counts = colSums(m), detected = colSums(m > 0), zero = colSums(m == 0))
}) |> left_join(cells, by = "cell")
stopifnot(all(per_cell$detected + per_cell$zero == per_cell$features_in_table))

# ---- Per cell line -------------------------------------------------------------------------------
features_by_line <- imap_dfr(libraries, function(m, lib) {
  map_dfr(PSCSR_CELL_LINES, function(line) {
    x <- m[, intersect(cells$cell[cells$cell_line == line], colnames(m)), drop = FALSE] > 0
    tibble(library = lib, cell_line = line,
           never_detected_in_line = sum(rowSums(x) == 0),
           detected_in_every_cell = sum(rowSums(x) == ncol(x)))
  })
})
per_line <- per_cell |>
  group_by(library, cell_line) |>
  summarise(cells = n(), features_in_table = first(features_in_table),
            counts_median = median(total_counts), counts_mean = mean(total_counts),
            counts_q25 = quantile(total_counts, 0.25), counts_q75 = quantile(total_counts, 0.75),
            detected_median = median(detected), detected_mean = mean(detected),
            detected_q25 = quantile(detected, 0.25), detected_q75 = quantile(detected, 0.75),
            zero_median = median(zero), zero_mean = mean(zero), .groups = "drop") |>
  left_join(features_by_line, by = c("library", "cell_line")) |>
  mutate(cell_line = factor(cell_line, PSCSR_CELL_LINES)) |> arrange(library, cell_line)

write_tsv(per_cell, file.path(out_dir, "qc_per_cell.tsv"))
write_tsv(per_line, file.path(out_dir, "qc_per_cell_line.tsv"))
for (lib in names(libraries)) {
  cat(sprintf("\n=== %s library (%d features in the table) ===\n", lib, per_line$features_in_table[per_line$library == lib][1]))
  print(per_line |> filter(library == lib) |>
          transmute(cell_line, cells,
                    `counts median [IQR]`   = sprintf("%.0f [%.0f-%.0f]", counts_median, counts_q25, counts_q75),
                    `detected median [IQR]` = sprintf("%.0f [%.0f-%.0f]", detected_median, detected_q25, detected_q75),
                    `detected mean` = round(detected_mean), `zero median` = round(zero_median),
                    `never detected in line` = never_detected_in_line, `detected in every cell` = detected_in_every_cell),
        width = Inf)
}

# ---- Plots ---------------------------------------------------------------------------------------------
long <- per_cell |> mutate(cell_line = factor(cell_line, PSCSR_CELL_LINES)) |>
  pivot_longer(c(total_counts, detected, zero), names_to = "metric") |>
  mutate(metric = factor(metric, c("total_counts", "detected", "zero"),
                         c("Total counts per cell", "Features detected per cell", "Features with zero counts per cell")))
p <- ggplot(long, aes(cell_line, value, fill = cell_line)) +
  geom_violin(scale = "width", colour = NA, alpha = 0.7) +
  geom_boxplot(width = 0.15, outlier.size = 0.3, fill = "white") +
  facet_wrap(library ~ metric, scales = "free_y", ncol = 3) +
  scale_y_continuous(labels = scales::comma) +
  theme_bw(base_size = 11) + theme(legend.position = "none") +
  labs(x = NULL, y = NULL, title = "PSCSR-seq V2: per-cell QC by cell line")
ggsave(file.path(out_dir, "qc_per_cell_line.pdf"), p, width = 12, height = 7)
cat(sprintf("\nWrote qc_per_cell.tsv, qc_per_cell_line.tsv, qc_per_cell_line.pdf to %s\n", out_dir))
