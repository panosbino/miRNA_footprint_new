library(tidyverse)
setwd("~/Desktop/Projects/miRNA_footprint_new/")
df <- read.delim("./datasets/sc_miRNA_seq/mirna_count_table_v2_2nd_run.tsv")
colnames(df) <- df |> colnames() |> str_split_i(i=1, pattern = "_")

mir124 <- df[df$mirna == ">Hsa-Mir-124-P1-v1_3p", ] |>
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
             as.numeric()) + 6
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
df_norm <- df |> column_to_rownames("mirna")  |> edgeR::cpm() |> as.data.frame()

colnames(df_norm) <- df_norm |> colnames() |> str_split_i(i=1, pattern = "_")

mir124_norm <- df_norm[">Hsa-Mir-124-P1-v1_3p", ] |>
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
