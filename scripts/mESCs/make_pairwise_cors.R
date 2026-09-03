library(tidyverse)

setwd("~/Desktop/microIMP_for_Vi/mESC_inhouse")

# Load the Targets, counts and phenotipic info
targets_path = "./Data/Target_files/Targets_new.rds"
counts_path = "./Data/final_counts.rds"
pheno = read.delim("./Data/phenotype.csv", sep = ",") ; pheno$Exp[pheno$Exp == "WT"] = "Control" ; kos = pheno[pheno$Exp == "KO",]$cells ; controls = pheno[pheno$Exp == "Control",]$cells 
counts = readRDS(counts_path) ; counts = counts[,colnames(counts) %in% pheno$cells]
targets = readRDS(targets_path)

# Subset the counts to only include targets
counts_targets_only = counts[rownames(counts) %in% targets$ensembl_gene_id,]
# Subset the counts to only include control cells we don't need the KO. 
#But we could do one with the KO to show the lack of specific the network, maybe this is already in Marcel's paper?
counts_targets_only_control = counts_targets_only[,controls]
counts_targets_only_control = counts_targets_only_control %>% t() %>% as.data.frame()

# Write the combinations
a = expand.grid(rownames(counts_targets_only),rownames(counts_targets_only))
a$Var1 = a$Var1 %>% as.character()
a$Var2 = a$Var2 %>% as.character()

#Calc the correlations. Somewhere I found that bind_rows() is faster than the purrr:reduce(rbind)
corrrrs = lapply(1:nrow(a), function(i) {
  cor = cor.test(counts_targets_only_control[,a[i,1]],counts_targets_only_control[,a[i,2]], method = "spearman")
  return(c(a[i,1],a[i,2],cor$estimate,cor$p.value))
  print(i)
}) #%>% bind_rows()

#saveRDS(corrrrs,"~/Desktop/microIMP/mESC/Target_comp/correlation_list.rds")

system.time({cor_df = corrrrs %>% as.data.frame() %>% t() %>% as.data.frame()})

# Calculate p adjusted values
colnames(cor_df) = c("Gene1", "Gene2", "rho", "pval")
cor_df$padj = p.adjust(cor_df$pval,method = "BH")
rownames(cor_df) = seq(1:nrow(cor_df))
#cor_df = readRDS("~/Desktop/microIMP/mESC/Target_comp/correlation_df.rds")
cor_df$rho = cor_df$rho %>% as.numeric()
cor_df$pval = cor_df$pval %>% as.numeric()
cor_df$padj = cor_df$padj %>% as.numeric()
cor_df = cor_df %>% filter(rho < X)
#saveRDS(cor_df,"~/Desktop/microIMP/mESC/Target_comp/correlation_df.rds")

#ggplot() + geom_density(data = cor_df, mapping = aes(x = V3)) + theme_bw()









