library(tidyverse)
setwd("~/Desktop/Projects/miRNA_footprint_new/")
df <- read.delim("./datasets/sc_miRNA_seq/counts_clean.tsv")
colnames(df) <- df |> colnames() |> str_split_i(i=1, pattern = "_")
df <- df |> dplyr::select(-E04)  # failed well (all reads from one miRNA); excluded from all analyses

mir124 <- df[df$feature == "Hsa-Mir-124-P1-v1_3p", ] |>
  as.data.frame() |>
  t() |>
  as.data.frame()

mir124 <- mir124[-1, , drop = FALSE]
colnames(mir124) <- "value"

mir124 <- mir124 |>
  mutate(value = as.numeric(value)) |>
  arrange(desc(value)) |>
  rownames_to_column("wells") |>
  mutate(
    row = gsub("[0-9]", "", wells),
    col = (gsub("[A-Za-z]", "", wells) |>
             str_remove("0") |>
             as.numeric()) + 6,
    wells = paste0(row, col)  # physical well ID (columns 7-12), so it cannot collide with the zero wells (1-5)
  )

mir124_sorted <- mir124[order(mir124$row, mir124$col), ]

# Create the ordered factor levels based on the sorted order
row_levels <- mir124_sorted$row
mir124_sorted$row <- factor(mir124_sorted$row, levels = rev(row_levels |> unique()), ordered = TRUE)

mir124_sorted$col <- mir124_sorted$col 
col_levels <- mir124_sorted$col
mir124_sorted$col <- factor(mir124_sorted$col, levels = col_levels |> unique(), ordered = TRUE)

#mir124_sorted[mir124_sorted$wells == "D04",2] <- 40
plate_2nd_part <- ggplot() +
    geom_text(data = mir124_sorted, mapping = aes(x = col, y = row, label = value), size = 8) +
    theme_bw() +
    theme(axis.text = element_text(size = 20))

ggsave(plot = plate_2nd_part, filename = "./analysis/sc_miRNA_seq/plate_2nd_part.pdf", device = "pdf", width = 6, height = 4)

write.csv(mir124_sorted, "./analysis/sc_miRNA_seq/mir124_counts_v2_2nd.csv")

# normalize
library(edgeR)
df_norm <- df |> column_to_rownames("feature")  |> edgeR::cpm() |> as.data.frame()

# TMM normalization (edgeR): CPM using TMM-adjusted library sizes
dge <- df |> column_to_rownames("feature") |> as.matrix() |> DGEList()
dge <- calcNormFactors(dge, method = "TMM")
df_tmm <- cpm(dge) |> as.data.frame()

# DESeq2 normalization: counts divided by size factors.
# type = "poscounts" because no miRNA is non-zero in all wells, so the default
# median-of-ratios cannot compute its geometric means and stops with an error.
library(DESeq2)
count_mat <- df |> column_to_rownames("feature") |> as.matrix()
dds <- DESeqDataSetFromMatrix(countData = count_mat,
                              colData = data.frame(row.names = colnames(count_mat)),
                              design = ~ 1)
dds <- estimateSizeFactors(dds, type = "poscounts")
df_deseq2 <- counts(dds, normalized = TRUE) |> as.data.frame()

write.csv(df_tmm, "./datasets/sc_miRNA_seq/mirna_counts_TMM_v2_2nd.csv")
write.csv(df_deseq2, "./datasets/sc_miRNA_seq/mirna_counts_DESeq2_v2_2nd.csv")

# Full plate for TMM and DESeq2: sequenced wells (physical columns 7-12) plus the
# GFP-negative wells (physical columns 1-5), which were not sequenced and are set to 0
mir124_full_plate <- function(norm_df) {
  sequenced <- tibble(wells = colnames(norm_df),
                      value = as.numeric(norm_df["Hsa-Mir-124-P1-v1_3p", ])) |>
    mutate(row = str_sub(wells, 1, 1),
           col = as.integer(str_sub(wells, 2, 3)) + 6,
           wells = paste0(row, col))
  zeros <- tibble(row = rep(LETTERS[1:8], each = 5),
                  col = rep(1:5, times = 8),
                  value = 0) |>
    mutate(wells = paste0(row, col))
  bind_rows(zeros, sequenced) |>
    dplyr::select(wells, value, row, col) |>
    arrange(row, col)
}

write.csv(mir124_full_plate(df_tmm), "./datasets/sc_miRNA_seq/mir124_counts_TMM_full.csv")
write.csv(mir124_full_plate(df_deseq2), "./datasets/sc_miRNA_seq/mir124_counts_DESeq2_full.csv")

colnames(df_norm) <- df_norm |> colnames() |> str_split_i(i=1, pattern = "_")

mir124_norm <- df_norm["Hsa-Mir-124-P1-v1_3p", ] |>
  as.data.frame() |>
  t() |>
  as.data.frame()

colnames(mir124_norm) <- "value"

mir124_norm <- mir124_norm |>
  mutate(value = as.numeric(value)) |>
  arrange(desc(value)) |>
  rownames_to_column("wells") |>
  mutate(
    row = gsub("[0-9]", "", wells),
    col = (gsub("[A-Za-z]", "", wells) |>
             str_remove("0") |>
             as.numeric()) + 6,
    value = value |> round(digits = 0),
    wells = paste0(row,col)
  )

mir124_norm_sorted <- mir124_norm[order(mir124_norm$row, mir124_norm$col), ]

# Create the ordered factor levels based on the sorted order
row_levels <- mir124_norm_sorted$row
mir124_norm_sorted$row <- factor(mir124_norm_sorted$row, levels = rev(row_levels |> unique()), ordered = TRUE)

mir124_norm_sorted$col <- mir124_norm_sorted$col 
col_levels <- mir124_norm_sorted$col
mir124_norm_sorted$col <- factor(mir124_norm_sorted$col, levels = col_levels |> unique(), ordered = TRUE)

#mir124_norm_sorted[mir124_norm_sorted$wells == "D04",2] <- 40
plate_2nd_part <- ggplot() +
  geom_text(data = mir124_norm_sorted, mapping = aes(x = col, y = row, label = value), size = 8) +
  theme_bw() +
  theme(axis.text = element_text(size = 20))

plate_2nd_part
mir124_norm_sorted

ggsave(plot = plate_2nd_part, filename = "./analysis/sc_miRNA_seq/plate_2nd_part.pdf", device = "pdf", width = 6, height = 4)

write.csv(mir124_norm_sorted, "./analysis/sc_miRNA_seq/mir124_norm_counts_v2_2nd.csv")

wells_zero <- paste0(
  rep(LETTERS[1:8], each = 5),
  sprintf("%02d", rep(1:5, times = 8))
)

df_zeros <- data.frame(
  "value" = rep(0, length(wells_zero)),
  row.names = wells_zero,
  check.names = FALSE
)

df_zeros <- df_zeros |> rownames_to_column("wells") |> mutate(row = wells |> str_sub(start = 1, end = 1),
                                                              col = wells |> str_sub(start = 3, end = 3))
mir124_norm_sorted <- rbind(df_zeros,mir124_norm_sorted)

mir124_norm_sorted <- mir124_norm_sorted |>
  mutate(
    wells = if_else(
      substr(wells, 2, 2) == "0",
      paste0(substr(wells, 1, 1), substr(wells, 3, nchar(wells))),
      wells
    )
  )


mir124_sorted <- rbind(df_zeros,mir124_sorted)

mir124_sorted <- mir124_sorted |>
  mutate(
    wells = if_else(
      substr(wells, 2, 2) == "0",
      paste0(substr(wells, 1, 1), substr(wells, 3, nchar(wells))),
      wells
    )
  )

write.csv(mir124_sorted, "./datasets/sc_miRNA_seq/mir124_counts_full.csv")
write.csv(mir124_norm_sorted, "./datasets/sc_miRNA_seq/mir124_counts_normalized_full.csv")
