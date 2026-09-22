library(tidyverse)


setwd("~/Desktop/Projects/miRNA_footprint_new/")
sc_norm <- read.delim("analysis/sc_miRNA_seq/mir124_counts_normalized_full.csv", sep = ",", row.names = 1)

gfp <- read_delim("datasets/HEK_SS3/FACS_data/Plate4/Plate4_GFP_fluorescence.tsv")
colnames(gfp) <- c("wells","gfp","row","col")

our_cor_df <- merge(sc_norm, gfp, by = c("wells","row","col"))
#our_cor_df <- our_cor_df[our_cor_df$eGFP < gfp_upper_cutoff, ]
our_cor_df <- our_cor_df[is.finite(our_cor_df$value) & is.finite(our_cor_df$gfp), ]
our_cor_df <- our_cor_df |> mutate(gfp = case_when(gfp <= 0 ~ 0, .default = gfp))
our_cor <- cor.test(our_cor_df$value, our_cor_df$gfp, method = "spearman", exact = FALSE)
cat(sprintf("Our method (scran, TargetScan): rho = %.3f, n = %d\n", our_cor$estimate, nrow(our_cor_df)))

library(ggplot2)
library(patchwork)

# ---------------------------------------------------------
# 1. Data
# ---------------------------------------------------------

# Original data: includes cells where value = 0
df_all <- our_cor_df

# Remove cells where miR-124 value = 0
df_nozero <- our_cor_df %>%
  dplyr::filter(value != 0)


# ---------------------------------------------------------
# 2. Calculate Spearman correlations
# ---------------------------------------------------------

# --- All cells ---

cor_all_1 <- cor(
  log10(df_all$gfp + 1),
  log10(df_all$value + 1),
  method = "spearman",
  use = "complete.obs"
)

cor_all_2 <- cor(
  log10(df_all$gfp + 1),
  df_all$value,
  method = "spearman",
  use = "complete.obs"
)

cor_all_3 <- cor(
  df_all$gfp,
  df_all$value,
  method = "spearman",
  use = "complete.obs"
)


# --- Cells with value != 0 ---

cor_nozero_1 <- cor(
  log10(df_nozero$gfp + 1),
  log10(df_nozero$value + 1),
  method = "spearman",
  use = "complete.obs"
)

cor_nozero_2 <- cor(
  log10(df_nozero$gfp + 1),
  df_nozero$value,
  method = "spearman",
  use = "complete.obs"
)

cor_nozero_3 <- cor(
  df_nozero$gfp,
  df_nozero$value,
  method = "spearman",
  use = "complete.obs"
)


# ---------------------------------------------------------
# 3. Plots - ALL CELLS
# ---------------------------------------------------------

p1_all <- ggplot(
  df_all,
  aes(x = log10(gfp + 1), y = log10(value + 1))
) +
  geom_point(fill = "lightblue", size = 3, shape = 21, color = "black") +
  labs(
    title = paste0(
      "Log10 GFP vs. Log10 miR-124 counts\n",
      "Spearman \u03c1 = ", round(cor_all_1, 3)
    ),
    x = "log10 GFP",
    y = "log10 miR-124 normalized counts"
  ) +
  theme_bw(base_size = 13) +
  theme(
    plot.title = element_text(hjust = 0.5)
  )


p2_all <- ggplot(
  df_all,
  aes(x = log10(gfp + 1), y = value)
) +
  geom_point(fill = "hotpink", size = 3, shape = 21, color = "black") +
  labs(
    title = paste0(
      "Log10 GFP vs. miR-124 counts\n",
      "Spearman \u03c1 = ", round(cor_all_2, 3)
    ),
    x = "log10 GFP",
    y = "miR-124 normalized counts"
  ) +
  theme_bw(base_size = 13) +
  theme(
    plot.title = element_text(hjust = 0.5)
  )


p3_all <- ggplot(
  df_all,
  aes(x = gfp, y = value)
) +
  geom_point(fill = "lightgreen", size = 3, shape = 21, color = "black") +
  labs(
    title = paste0(
      "GFP vs. miR-124 counts\n",
      "Spearman \u03c1 = ", round(cor_all_3, 3)
    ),
    x = "GFP",
    y = "miR-124 normalized counts"
  ) +
  theme_bw(base_size = 13) +
  theme(
    plot.title = element_text(hjust = 0.5)
  )


# ---------------------------------------------------------
# 4. Plots - value != 0
# ---------------------------------------------------------

p1_nozero <- ggplot(
  df_nozero,
  aes(x = log10(gfp + 1), y = log10(value + 1))
) +
  geom_point(fill = "lightblue", shape = 24, size = 3, color = "black") +
  labs(
    title = paste0(
      "Log10 GFP vs. Log10 miR-124 counts\n",
      "value \u2260 0 | Spearman \u03c1 = ", round(cor_nozero_1, 3)
    ),
    x = "log10 GFP",
    y = "log10 miR-124 normalized counts"
  ) +
  theme_bw(base_size = 13) +
  theme(
    plot.title = element_text(hjust = 0.5)
  )


p2_nozero <- ggplot(
  df_nozero,
  aes(x = log10(gfp + 1), y = value)
) +
  geom_point(fill = "hotpink", shape = 24, size = 3, color = "black") +
  labs(
    title = paste0(
      "Log10 GFP vs. miR-124 counts\n",
      "value \u2260 0 | Spearman \u03c1 = ", round(cor_nozero_2, 3)
    ),
    x = "log10 GFP",
    y = "miR-124 normalized counts"
  ) +
  theme_bw(base_size = 13) +
  theme(
    plot.title = element_text(hjust = 0.5)
  )


p3_nozero <- ggplot(
  df_nozero,
  aes(x = gfp, y = value)
) +
  geom_point(fill = "lightgreen", shape = 24, size = 3, color = "black") +
  labs(
    title = paste0(
      "GFP vs. miR-124 counts\n",
      "value \u2260 0 | Spearman \u03c1 = ", round(cor_nozero_3, 3)
    ),
    x = "GFP",
    y = "miR-124 normalized counts"
  ) +
  theme_bw(base_size = 13) +
  theme(
    plot.title = element_text(hjust = 0.5)
  )


# ---------------------------------------------------------
# 5. Combine all 6 plots
# ---------------------------------------------------------

combined_plot <- (
  p1_all + p2_all + p3_all +
    p1_nozero + p2_nozero + p3_nozero
) +
  plot_layout(ncol = 3)

combined_plot

ggsave(plot = combined_plot, width = 12, height = 8,path = "./analysis/sc_miRNA_seq/", filename = "sc_seq_GFP_cor.png", device = "png")
ggsave(plot = combined_plot, width = 12, height = 8,path = "./analysis/sc_miRNA_seq/", filename = "sc_seq_GFP_cor.pdf", device = "pdf")

