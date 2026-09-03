library(tidyverse)
library(igraph)
library(ggraph)
library(tidygraph)

# system.time({ test = bind_rows(merged_list_short) })
# system.time({ test = merged_list_short %>% purrr::reduce(rbind) })
# system.time({ test = rbind.data.frame() })

#Load all the date
# this is the pairwise corralation df with the pvals and rhos
cor_df = readRDS("~/Desktop/microIMP/mESC/Target_comp/correlation_df.rds")
cor_df$padj_bonf = p.adjust(cor_df$pval,method = "bonferroni")
pheno = read.delim("~/Desktop/microIMP/mESC/Figure1/data/phenotype.csv", sep = ",") ; pheno$Exp[pheno$Exp == "WT"] = "Control" ; kos = pheno[pheno$Exp == "KO",]$cells ; controls = pheno[pheno$Exp == "Control",]$cells 
counts = readRDS("~/Desktop/microIMP/mESC/Figure1/data/final_counts.rds") ; counts = counts[,colnames(counts) %in% pheno$cells]
targets = readRDS("~/Desktop/microIMP/mESC/Figure1/data/Targets/Targets_new.rds")
targets_short = targets %>% dplyr::select(Cumulative.weighted.context...score, ensembl_gene_id)
targets_short = targets_short %>% group_by(ensembl_gene_id) %>% summarise_all(min)
targets_short = targets_short %>% arrange(Cumulative.weighted.context...score) %>% head(n = 1000)

calculate_activity = function(counts_df, targets, mESC_data = T) {
  cell_sums <- counts_df[rownames(counts_df) %in% targets$ensembl_gene_id,] %>% colSums() %>% as.data.frame()
  if (mESC_data) {cell_sums_norm = cell_sums/mean(cell_sums[kos,])
  } else {
    cell_sums_norm = cell_sums/mean(cell_sums$.)
  } 
  cell_sums_norm$. = ifelse(cell_sums_norm$. == 0, cell_sums_norm %>% filter(. > 0) %>% min(),cell_sums_norm$. )
  activity <- -log2(cell_sums_norm) ; activity <- merge(activity, pheno, by.x = 0, by.y = "cells")
  ko_activ = activity[activity$Exp == "KO", ]; control_activ = activity[activity$Exp == "Control", ]
  d = (mean(ko_activ$.) - mean(control_activ$.))/sd(ko_activ$.)
  overlap = (2*pnorm(-abs(d)/2))*100
  diff <- (activity[activity$Exp == "Control",]$. %>% mean()) - (activity[activity$Exp == "KO",]$. %>% mean())
  return(c(activity,overlap,diff))
}


# changed the 0 to the min pval to be able to do calc, might not be necessary
cor_df$padj = ifelse(cor_df$padj == 0,  10^-120, cor_df$padj )

# Filter the table 
cor_df_1000 = cor_df[cor_df$Gene1 %in% targets_short$ensembl_gene_id,]
cor_df_1000 = cor_df_1000[cor_df_1000$Gene2 %in% targets_short$ensembl_gene_id,]
cor_df_filtered = cor_df_1000 %>% filter(padj_bonf < 0.05 & rho > 0)

#ggplot() + geom_point(data = cor_df %>% filter(rho > 0), mapping =  aes(x = padj, y = rho), size = 0.5) + theme_bw() + scale_x_log10()


# cor_df_filt = cor_df %>% filter(abs(rho) > 0.1)
# cor_df_filt = cor_df_filt[cor_df_filt$Gene1 %in% targets_short$ensembl_gene_id,]
# cor_df_filt = cor_df_filt[cor_df_filt$Gene2 %in% targets_short$ensembl_gene_id,]


min(cor_df$rho)
cor_df$rho = cor_df$rho %>% as.numeric()
cor_df$cor_sign = ifelse(cor_df$rho < 0, "Negative", "Positive")

ggplot() + geom_boxplot(data = cor_df, aes(x = cor_sign,y = abs(rho))) + theme_bw() + ylim(0,1)

rhos = seq(min(cor_df$rho), max(cor_df$rho), 0.001)
res_rho = lapply(rhos, function(i){
  #Create the network object
  if (i > 0){
    cor_df_filtered = cor_df %>% filter(rho > i) 
  } else {
    cor_df_filtered = cor_df%>% filter(rho < i)
  }
  network = graph_from_data_frame(d=cor_df_filtered, directed=F) %>% as_tbl_graph()

  #Cluster
  a = network %>% activate(nodes) %>% mutate(community = as.factor(group_components())) %>% as.data.frame()
  covaring_genes = a %>% filter(community == 1) %>% select(-community)
  #covaring_genes = unique(c(unique(cor_df_filtered$Gene1),unique(cor_df_filtered$Gene1))) %>% as.data.frame()
  #covaring_genes$dummy = "dummy"
  colnames(covaring_genes)[1] = "ensembl_gene_id"
  activity = calculate_activity(counts_df = counts, targets = covaring_genes)
  activity_matrix = activity[1:3] %>% as.data.frame(); rownames(activity_matrix) = activity_matrix$Row.names; activity_matrix = activity_matrix[,2:3]; colnames(activity_matrix) = c("Activity", "Pheno")
  overlap = activity[[4]]
  diff = activity[[5]]
  print(i)
  return(c(i,overlap,diff,nrow(covaring_genes)))
}) %>% purrr::reduce(rbind) %>% as.data.frame()

colnames(res_rho) = c("rho", "overlap %", "distance","N_targets")
rownames(res_rho) = seq(1:nrow(res_rho))

ggplot() + geom_point(data = res_rho, mapping =  aes(x = rho, y = `overlap %`), size = 0.5) + theme_bw()# + ylim(0,10) #+ scale_x_log10() + 
ggplot() + geom_point(data = res_rho, mapping =  aes(x = rho, y = `distance`), size = 0.5) + theme_bw() #+ scale_x_log10()
ggplot() + geom_point(data = res_rho, mapping =  aes(x = rho, y = `distance`), size = 0.5) + theme_bw() #+ scale_x_log10()



tmp = cor_df %>% filter(rho > 0.5) 
tmp2 = c(tmp$Gene1,tmp)
targets_05 = targets[targets$ensembl_gene_id %in% covaring_genes$ensembl_gene_id,]
targets_048 = targets[targets$ensembl_gene_id %in% covaring_genes$ensembl_gene_id,]


ggplot() + geom_density(data = cor_df,aes(x = rho)) + theme_bw()

## From here I tried different things that are not well documented.

# Try the graph approach

ggraph(network, layout = 'fr') + 
  #geom_edge_link() + 
  geom_node_point(size = 1) + 
  theme(legend.position = 'bottom')
intersect(a$name,b$name) %>% length()


c = b[!b$name %in% a$name,]

i =0.54
if (i > 0){
  cor_df_filtered = cor_df_filt %>% filter(rho > i) 
} else {
  cor_df_filtered = cor_df_filt %>% filter(rho < i)
}
network1 = graph_from_data_frame(d=cor_df_filtered, directed=F) %>% as_tbl_graph()

#Cluster
a1 = network1 %>% activate(nodes) %>% mutate(community = as.factor(group_components())) %>% as.data.frame()


i =-0.54
if (i > 0){
  cor_df_filtered = cor_df_filt %>% filter(rho > i) 
} else {
  cor_df_filtered = cor_df_filt %>% filter(rho < i)
}
network2 = graph_from_data_frame(d=cor_df_filtered, directed=F) %>% as_tbl_graph()

#Cluster
a2 = network2 %>% activate(nodes) %>% mutate(community = as.factor(group_components())) %>% as.data.frame()

a3 = a1[!a1$name %in% a2$name,]
a4 = a1[a1$name %in% a2$name,]

colnames(covaring_genes)[1] = "ensembl_gene_id"
covaring_genes = a3$name %>% as.data.frame()
covaring_genes$dummy = "dummy"
activity = calculate_activity(counts_df = counts, targets = covaring_genes)
activity_matrix = activity[1:3] %>% as.data.frame(); rownames(activity_matrix) = activity_matrix$Row.names; activity_matrix = activity_matrix[,2:3]; colnames(activity_matrix) = c("Activity", "Pheno")
overlap = activity[[4]]
diff = activity[[5]]
write(paste0(a4$name,","), "~/Desktop/go.txt")


targets_co = targets[targets$ensembl_gene_id %in% a4$name,]


ggplot() + geom_point()



res_pval = lapply(c(5 %o% 10^(-1:-140)), function(i){
  #Create the network object
  cor_df_filtered = cor_df_filt %>% filter(padj_bonf < i & rho > 0) 
  network = graph_from_data_frame(d=cor_df_filtered, directed=F) %>% as_tbl_graph()
  #Cluster
  a = network %>% activate(nodes) %>% mutate(community = as.factor(group_components())) %>% as.data.frame()
  covaring_genes = a %>% filter(community == 1) %>% select(-community)
  covaring_genes$dummy = "dummy"
  colnames(covaring_genes)[1] = "ensembl_gene_id"
  activity = calculate_activity(counts_df = counts, targets = covaring_genes)
  activity_matrix = activity[1:3] %>% as.data.frame(); rownames(activity_matrix) = activity_matrix$Row.names; activity_matrix = activity_matrix[,2:3]; colnames(activity_matrix) = c("Activity", "Pheno")
  overlap = activity[[4]]
  diff = activity[[5]]
  return(c(i,overlap,diff,nrow(covaring_genes)))
}) %>% purrr::reduce(rbind) %>% as.data.frame()

colnames(res_pval) = c("padj", "overlap %", "distance")
rownames(res_pval) = seq(1:nrow(res_pval))
res_pval$source = "pval_changing"

ggplot() + geom_density(data = cor_df %>% filter(padj_bonf < 0.05), aes(x = log10(padj_bonf))) + theme_bw()

ggplot() + geom_point(data = res_pval, mapping =  aes(x = padj, y = `overlap %`), size = 0.5) + theme_bw() + scale_x_log10() + ylim(0,10)
ggplot() + geom_point(data = res_pval, mapping =  aes(x = padj, y = `distance`), size = 0.5) + theme_bw() + scale_x_log10()
ggplot() + geom_point(data = cor_df %>% filter(rho > 0), mapping =  aes(x = padj_bonf, y = rho), size = 0.5) + theme_bw() + scale_x_log10()



##########################


tmp = merge(x = cor_df, y = targets_short, by.x = "Gene1", by.y = "ensembl_gene_id", all.x = T)
tmp2 = merge(x = tmp, y = targets_short, by.x = "Gene2", by.y = "ensembl_gene_id", all.x = T)
tmp2$score_sum = (tmp2$Cumulative.weighted.context...score.x + tmp2$Cumulative.weighted.context...score.y)
tmp2$density = get_density(x = tmp2$score_sum, y = tmp2$rho)


ggplot() + geom_point(data = tmp2, mapping =  aes(x = score_sum, y = rho), size = 0.5) + theme_bw()# + scale_x_log10() + ylim(0,10)

#####################3

p_cov <- ggplot(activity_matrix, aes(x = Activity, fill = Pheno)) +
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
  xlim(c(-1.4,1.4)) +
  #xlim(c(min(activity_matrix$Activity)-0.2,max(activity_matrix$Activity)+0.2)) +
  labs(title = paste0( "Top 200     Overlap = ", round(overlap, digits = 3), "%   ", "Mean difference = ", round(diff, digits = 2)), x="Activity ", y = "Density")
p_cov


