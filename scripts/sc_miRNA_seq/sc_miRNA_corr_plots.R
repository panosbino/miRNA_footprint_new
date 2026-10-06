library(tidyverse)
library(patchwork)

setwd("~/Desktop/Projects/miRNA_footprint_new/")

# Full-plate miR-124 files written by prepare_miRNA_sc_mir124_counts.R
norm_files <- c(
  CPM    = "datasets/sc_miRNA_seq/mir124_counts_normalized_full.csv",
  TMM    = "datasets/sc_miRNA_seq/mir124_counts_TMM_full.csv",
  DESeq2 = "datasets/sc_miRNA_seq/mir124_counts_DESeq2_full.csv"
)

# Induced cells = sequenced wells, physical columns 7-12 (defined by plate layout,
# not by GFP). Columns 1-5 are uninduced wells with synthetic zeros: excluded.
INDUCED_COLS <- 7:12

# Detected cells: at least MIN_RAW_READS raw miR-124 reads. Removes the zero cells and
# the 1-read cell (B8); the next-lowest cell has 13 reads, so any value 2-13 gives the same set.
MIN_RAW_READS <- 2
raw_reads <- read.delim("datasets/sc_miRNA_seq/mir124_counts_full.csv", sep = ",", row.names = 1) |>
  dplyr::select(wells, raw_reads = value)

gfp <- read_delim("datasets/HEK_SS3/FACS_data/Plate4/Plate4_GFP_fluorescence.tsv")
colnames(gfp) <- c("wells", "gfp", "row", "col")

# Spearman and Pearson correlation with p-values
correlate <- function(x, y) {
  sp <- cor.test(x, y, method = "spearman", exact = FALSE)
  pe <- cor.test(x, y, method = "pearson")
  tibble(spearman_rho = unname(sp$estimate), spearman_p = sp$p.value,
         pearson_r = unname(pe$estimate), pearson_p = pe$p.value,
         n = sum(complete.cases(x, y)))
}

# One scatter plot in the style of the original script
plot_cor <- function(df, x, y, title, xlab, ylab, fill, shape, color, cor_res, prefix = "") {
  ggplot(df, aes(x = {{ x }}, y = {{ y }})) +
    geom_point(fill = fill, size = 3, shape = shape, color = color) +
    labs(
      title = paste0(title, "\n", prefix, "n = ", cor_res$n, "\n",
                     "Spearman \u03c1 = ", round(cor_res$spearman_rho, 3), ", p = ", signif(cor_res$spearman_p, 2), "\n",
                     "Pearson r = ", round(cor_res$pearson_r, 3), ", p = ", signif(cor_res$pearson_p, 2)),
      x = xlab,
      y = ylab
    ) +
    theme_bw(base_size = 13) +
    theme(plot.title = element_text(hjust = 0.5, size = 11))
}

cor_summary <- list()

for (norm in names(norm_files)) {

  # ---------------------------------------------------------
  # 1. Data
  # ---------------------------------------------------------
  sc_norm <- read.delim(norm_files[[norm]], sep = ",", row.names = 1)

  our_cor_df <- merge(sc_norm, gfp, by = c("wells", "row", "col"))
  our_cor_df <- our_cor_df[our_cor_df$col %in% INDUCED_COLS, ]
  our_cor_df <- our_cor_df[is.finite(our_cor_df$value) & is.finite(our_cor_df$gfp), ]
  our_cor_df <- our_cor_df |> mutate(gfp = case_when(gfp <= 0 ~ 0, .default = gfp))

  # Original data: includes cells where value = 0
  df_all <- our_cor_df

  # Keep cells with miR-124 detected (>= MIN_RAW_READS raw reads)
  df_nozero <- our_cor_df |>
    left_join(raw_reads, by = "wells") |>
    dplyr::filter(raw_reads >= MIN_RAW_READS)

  # ---------------------------------------------------------
  # 2. Correlations: log10 GFP vs log10 miR-124
  # ---------------------------------------------------------
  cor_all_1      <- correlate(log10(df_all$gfp + 1),    log10(df_all$value + 1))
  cor_nozero_1   <- correlate(log10(df_nozero$gfp + 1), log10(df_nozero$value + 1))

  cor_summary[[norm]] <- bind_rows(
    cor_all_1      |> mutate(cells = "all",        comparison = "log10 GFP vs log10 miR-124"),
    cor_nozero_1   |> mutate(cells = paste0("reads >= ", MIN_RAW_READS), comparison = "log10 GFP vs log10 miR-124")
  ) |>
    mutate(normalization = norm, .before = 1)

  cat(sprintf("%s, induced cells: Spearman rho = %.3f (p = %.3g), Pearson r = %.3f (p = %.3g), n = %d\n",
              norm, cor_all_1$spearman_rho, cor_all_1$spearman_p,
              cor_all_1$pearson_r, cor_all_1$pearson_p, cor_all_1$n))

  # ---------------------------------------------------------
  # 3. Plot - ALL CELLS
  # ---------------------------------------------------------
  ylab_log <- paste0("log10 miR-124 (", norm, ")")

  p1_all <- plot_cor(df_all, log10(gfp + 1), log10(value + 1),
                     paste0("Log10 GFP vs. Log10 miR-124, ", norm), "log10 GFP", ylab_log,
                     "#7996ec", 21, "#406ae4", cor_all_1)

  # ---------------------------------------------------------
  # 4. Plot - miR-124 detected (>= MIN_RAW_READS reads)
  # ---------------------------------------------------------
  p1_nozero <- plot_cor(df_nozero, log10(gfp + 1), log10(value + 1),
                        paste0("Log10 GFP vs. Log10 miR-124, ", norm), "log10 GFP", ylab_log,
                        "lightblue", 24, "black", cor_nozero_1, prefix = paste0("miR-124 reads \u2265 ", MIN_RAW_READS, " | "))

  # ---------------------------------------------------------
  # ---------------------------------------------------------
  # 5. Combine the two log-log plots
  # ---------------------------------------------------------
  combined_plot <- (p1_all + p1_nozero) +
    plot_layout(ncol = 2) +
    plot_annotation(title = paste0("miR-124 vs GFP, induced cells, ", norm, " normalization"))

  print(combined_plot)

  ggsave(plot = combined_plot, width = 9, height = 5.5, path = "./analysis/sc_miRNA_seq/",
         filename = paste0("sc_seq_GFP_cor_induced_", norm, ".png"), device = "png")
  ggsave(plot = combined_plot, width = 9, height = 5.5, path = "./analysis/sc_miRNA_seq/",
         filename = paste0("sc_seq_GFP_cor_induced_", norm, ".pdf"), device = cairo_pdf)
  ggsave(plot = p1_all, width = 4.5, height = 5, path = "./analysis/sc_miRNA_seq/",
         filename = paste0("sc_seq_GFP_cor_log_log_induced_", norm, ".pdf"), device = cairo_pdf)
}

# ---------------------------------------------------------
# 6. Summary table of all correlations
# ---------------------------------------------------------
cor_summary <- bind_rows(cor_summary) |>
  mutate(across(c(spearman_rho, pearson_r), \(x) round(x, 3)),
         across(c(spearman_p, pearson_p), \(x) signif(x, 3)))
print(cor_summary, n = Inf)
write_csv(cor_summary, "./analysis/sc_miRNA_seq/sc_seq_GFP_cor_induced_summary.csv")
