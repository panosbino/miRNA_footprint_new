library(tidyverse)

calculate_activity = function(counts_df, targets, mESC_data = T) {
  cell_sums <- counts_df[rownames(counts_df) %in% targets$ensembl_gene_id,] %>% colSums() %>% as.data.frame()
  if (mESC_data) {cell_sums_norm = cell_sums/mean(cell_sums[kos,])
  } else {
    cell_sums_norm = cell_sums/mean(cell_sums$.)
  } 
  activity <- -log2(cell_sums_norm) ; activity <- merge(activity, pheno, by.x = 0, by.y = "cells")
  ko_activ = activity[activity$Exp == "KO", ]; control_activ = activity[activity$Exp == "Control", ]
  d = (mean(ko_activ$.) - mean(control_activ$.))/sd(ko_activ$.)
  overlap = (2*pnorm(-abs(d)/2))*100
  diff <- (activity[activity$Exp == "Control",]$. %>% mean()) - (activity[activity$Exp == "KO",]$. %>% mean())
  return(c(activity,overlap,diff))
}

setwd("~/Desktop/Projects/miRNA_footprint_new/")
pheno = read.delim("./datasets/mESCs/phenotype.csv", sep = ",")
# Change the name from WT to controls as it is more accurate
pheno$Exp[pheno$Exp == "WT"] = "Control" 
# Find the KO and control cells
kos = pheno[pheno$Exp == "KO",]$cells ; controls = pheno[pheno$Exp == "Control",]$cells 
# Read the count data
counts = readRDS("datasets/mESCs/final_counts.rds") ; counts = counts[,colnames(counts) %in% pheno$cells]
# Load the pre-processed targets
targets = readRDS("resources/Targets_combined_top3_families_mESCs.rds")
# Combine scores of genes that are targets of more than one miRNA
targets = targets %>% group_by(ensembl_gene_id) %>% mutate(score = sum(Cumulative.weighted.context...score)) %>% ungroup() %>% distinct(pick(ensembl_gene_id),.keep_all = T)

activity = calculate_activity(counts_df = counts, targets = targets, mESC_data = T)
activity_matrix = activity[1:3] %>% as.data.frame(); rownames(activity_matrix) = activity_matrix$Row.names; activity_matrix = activity_matrix[,2:3]; colnames(activity_matrix) = c("Activity", "Pheno")
overlap = activity[[4]]
diff = activity[[5]]

p1 <- ggplot(activity_matrix, aes(x = Activity, fill = Pheno)) +
  geom_density(color="black", alpha=0.8, outline.type = "both", size = 1.3) +
  theme(panel.border = element_rect(colour = "black", linewidth = 2, fill = NA),
        panel.background = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major.y = element_line(colour = "black", linetype = "dashed", linewidth = 0.2),
        panel.grid.major.x = element_blank(),
        axis.title.x = element_text( size = 16),
        axis.title.y = element_text( size = 16),
        plot.title = element_text(size = 18, hjust = 0.5),
        legend.title = element_blank())+
  scale_fill_manual(values = c("royalblue4","lightblue1")) +
  #xlim(c(-0.2,0.2)) +
  #xlim(c(min(activity_matrix$Activity)-0.2,max(activity_matrix$Activity)+0.2)) +
  labs(title = paste0( "All targets     Overlap = ", round(overlap, digits = 3), "%   ", "Mean difference = ", round(diff, digits = 2)), x="Activity ", y = "Density")
p1












