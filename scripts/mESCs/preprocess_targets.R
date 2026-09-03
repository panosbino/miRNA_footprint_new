# Newer versions of dplyr don't work with biomart currently
#devtools::install_version("dbplyr", version = "2.3.4")

library(dplyr)
library(stringr)
library(biomaRt)

setwd("~/Desktop/microIMP_for_Vi/mESC_inhouse")
ensembl_mm = useEnsembl(biomart="ensembl", dataset="mmusculus_gene_ensembl")
#library(tidyverse)
# First we need to combine the target files of the different miRNA families.
# Some files have different columns because of different conservation category i think. 
# to combine them we need to keep only the common colums

a = read.delim("./Data/Target_files/TargetScan8.0__miR-17-5p_20-5p_93-5p_106-5p.predicted_targets.txt")
b = read.delim("./Data/Target_files/TargetScan8.0__miR-291-3p_294-3p_295-3p_302-3p.predicted_targets.txt")
c = read.delim("./Data/Target_files/TargetScan8.0__miR-292a-3p_467a-5p.predicted_targets.txt")

common_cols = intersect(colnames(a), colnames(b)) %>% intersect(colnames(c))
a = a %>% dplyr::select(all_of(common_cols))
b = b %>% dplyr::select(all_of(common_cols))
c = c %>% dplyr::select(all_of(common_cols))

# This is the combined file
targets = rbind(a,b) %>% rbind(c)
rm(common_cols,a,b,c)

#write_delim(targets, "~/Desktop/microIMP/mESC/Targets__NC_new.txt", delim = "\t")


# Then we need to change the ensemble transcript IDs to ensembl gene IDs to be able to match it with the count matrix in downstream analysis
#Remove the version part from ensemble transcript ID so it can match the the biomart file
first.word = function(my.string){
  unlist(str_split(my.string,fixed(".")))[1]
}
#Function that changes the ensembl transcript IDs to ensembl gene IDs (from targetscan)
switch_names_targets = function(targets) {
  no_genes = dim(targets)[1]
  genes = getBM(attributes=c('ensembl_gene_id','ensembl_transcript_id'), 
                 filters = 'ensembl_transcript_id', 
                 values = targets$Representative.transcript %>% unique(), 
                 mart = ensembl_mm)
  targets2 = merge(targets, genes, by.x = "Representative.transcript", by.y = "ensembl_transcript_id", all.x = TRUE)
  cat(dim(genes)[1] ,"ENSEMBL transcript IDs successfully changed ENSEMBL gene IDs, out of", no_genes, "total genes" )
  return(targets2)
}



main = function(targets) {
  targets$Representative.transcript = sapply(targets$Representative.transcript, first.word)
  targets_new = switch_names_targets(targets = targets)
  targets_new = arrange(targets_new, Cumulative.weighted.context...score)
  targets_new = targets_new[!is.na(targets_new$ensembl_gene_id),] 
  return(targets_new)
}

#path_to_targetscan_file = "./Data/Target_files/Targets__NC_new.txt"
targets = read.delim(path_to_targetscan_file, sep = "\t")
targets_new = main(targets = targets)

# Extract the name 
name = path_to_targetscan_file %>% str_split_i(pattern = "__", i = 2) %>% str_split_i(pattern = fixed("."), i = 1)
# or just enter it if you don't load a file
name = "combined_NEW"
saveRDS(targets_new, paste0("./Data/Target_files/Targets_", name,".rds"))

