#Author: Panos
#Purpose: This script is used to preprocess the count matrix for the Rzepieal data
#The count matrix is loaded, cpm normalized, and the genes names are changed to ensemble gene IDs. 
#The targets are loaded, the transcript IDs are switched to gene IDs so they match the count matrix nnd the targets are arranged with respect to the target score. 
#Also the "good cellls" are saved, these are cells that have more than 47700 raw counts. This is based on the corresponding analysis.
#Finally, the normalized GFP counts are also saved.

library(Seurat)
library(tidyverse)
library(biomaRt)
library("ggpubr")
library(MASS)
library(viridis)

setwd("~/Desktop/microIMP_for_Vi/HEK_Rzepiela")

#Load the count matrix
data<- Read10X("./Data/filtered_feature_bc_matrix/")
data <- CreateSeuratObject(data)
#Get the matrix with the raw counts
raw_counts <- GetAssayData(data, slot = "counts") %>% as.matrix() %>% as.data.frame()
saveRDS(raw_counts, "raw_counts.rds")
counts_per_cell = raw_counts %>% colSums() %>% as.data.frame()
ggplot() + geom_histogram(data = counts_per_cell, mapping = aes(x = .), color = "black", fill = "#fbe5a0") + theme_bw()


# #Find the deeply sequenced cells, based on 47700 analysis
# good_cells <- raw_counts %>% colSums() 
# good_cells <- good_cells[good_cells > 47700] %>% names()
# saveRDS(good_cells, "./good_cells.rds")

#Normalize the count matrix, cpm
data <- NormalizeData(data, normalization.method = "RC", scale.factor = 1000000)
counts <- GetAssayData(data, slot = "data") %>% as.matrix() %>% as.data.frame()
#Get the GFP normalized counts
GFP_counts <- counts["eGFP",] %>% t() %>% as.data.frame()
#Make the biomart object 
ensembl_hs = useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl") ### Has to run before the function below
#Get the table with the ensembl gene IDs and the gene names 
genes <- getBM(attributes=c('ensembl_gene_id',"external_gene_name"), 
               filters ='external_gene_name', 
               values = rownames(counts), 
               mart = ensembl_hs)
#Change the names
counts <- merge(counts, genes, by.x = 0, by.y = "external_gene_name")
rownames(counts) <- counts[,dim(counts)[2]]
counts <- counts[,3:dim(counts)[2]-1]
saveRDS(counts, "./Data/counts.rds")
saveRDS(GFP_counts, "./Data/GFP_counts.rds")
