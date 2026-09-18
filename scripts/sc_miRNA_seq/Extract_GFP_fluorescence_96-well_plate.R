library(flowCore)
library(flowGraph)

setwd("~/Desktop/Projects/miRNA_footprint_new/")

plate <- "Plate4"

# 96-well plate dimensions 
plate_rows <- 8 
plate_columns <- 12 
# Numbers from 0 to 7 
nums <- 0:(plate_rows - 1) 
# Convert to letters: A = 0, B = 1, ..., H = 7 
letters_equiv <- LETTERS[nums + 1]

# Combine into a named vector or data frame
translation_table <- data.frame(
  Number = nums,
  Letter = letters_equiv
)

# Input path
path <- paste0("datasets/HEK_SS3/FACS_data/", plate, "/")

file_names <- system(
  paste0("ls ", path, "*.fcs"),
  intern = TRUE
)

dfs <- map_df(file_names, function(file_name) {
  x <- read.FCS(file_name, transformation = FALSE)
  
  df <- getIndexSort(x = x)
  
  df <- df |>
    inner_join(
      translation_table,
      by = join_by("XLoc" == "Number"),
      relationship = "many-to-one"
    ) |>
    dplyr::mutate(YLoc = YLoc + 1) |>
    dplyr::rename(
      Row = "Letter",
      Column = "YLoc",
      GFP = "FITC.A"
    ) |>
    dplyr::mutate(
      Well_ID = paste0(Row, Column)
    ) |>
    dplyr::select(
      Well_ID,
      GFP,
      Row,
      Column
    )
})

# Output paths
output_dir <- path

output_tsv <- file.path(
  output_dir,
  paste0(plate, "_GFP_fluorescence.tsv")
)

output_png <- file.path(
  output_dir,
  paste0(plate, "_GFP_fluorescence.png")
)

# Write fluorescence data
write_delim(
  x = dfs,
  file = output_tsv,
  delim = "\t"
)

# Plot
dfs$Row <- dfs$Row |>
  factor(levels = rev(unique(dfs$Row)))

plate_image <- dfs |>
  ggplot() +
  geom_point(
    aes(
      x = Column,
      y = Row,
      fill = log10(GFP)
    ),
    shape = 21,
    size = 6
  ) +
  theme_bw() +
  scale_fill_gradient2(
    low = "gray",
    mid = "white",
    high = "green",
    midpoint = 2
  ) +
  scale_x_continuous(
    breaks = seq(0, 24, 1)
  ) +
  theme(
    panel.grid.minor = element_blank()
  )

# Save plot
ggsave(
  plot = plate_image,
  path = output_dir,
  filename = paste0(plate, "_GFP_fluorescence.png"),
  device = "png",
  width = 8,
  height = 4
)

