library(flowCore)
library(flowGraph)

# Numbers from 0 to 15
nums <- 0:15

# Convert to letters: A = 0, B = 1, ..., P = 15
letters_equiv <- LETTERS[nums + 1]

# Combine into a named vector or data frame
translation_table <- data.frame(Number = nums, Letter = letters_equiv)

path <- "~/Desktop/Projects/miRNA_footprint_new/datasets/HEK_SS3/FACS_data/Plate5/"
file_names <- system(paste0("ls ",path,"*.fcs"), intern = T)

dfs <- map_df(file_names, function(file_name){
  x <- read.FCS(file_name, transformation=FALSE)
  df <- getIndexSort(x = x)
  df <- df |>
    inner_join(translation_table, by = join_by("XLoc" == "Number"), relationship = "many-to-one") |>
    dplyr::mutate(YLoc = YLoc + 1) |>
    dplyr::rename(Row = "Letter", Column = "YLoc", GFP = "FITC.A") |>
    dplyr::mutate(Well_ID = paste0(Row, Column)) |>
    dplyr::select(Well_ID, GFP,Row, Column) 
})

write_delim(x = dfs, file = "/Users/panagiotiskalogeropoulos/Desktop/Projects/microIMP/FACS_data/P105_GFP_fluorescence.tsv", delim = "\t")


dfs$Row <-dfs$Row |> factor(levels = c(rev(unique(dfs$Row))))
  
plate_image <- dfs |>
  ggplot() +
  geom_point(aes(x = Column, y = Row, fill = log10(GFP)), shape = 21, size = 6 ) +
  theme_bw() +
  scale_fill_gradient2(low = "gray", mid = "white", high = "green", midpoint = 2) +
  scale_x_continuous(breaks = seq(0, 24, 1)) +
  theme(panel.grid.minor = element_blank())


ggsave( plot = plate_image, path = "/Users/panagiotiskalogeropoulos/Desktop/Projects/microIMP/FACS_data/", filename = "P105_GFP_fluorescence.png", device = "png", width = 8, height = 4)


