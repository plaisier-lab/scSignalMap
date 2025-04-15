# docker run -it -v '/home/cplaisier/scSignalMap:/files' cplaisier/ccafv2_s5_extra

library(Seurat)
library(dplyr)

setwd('/files/test')

seurat_obj = readRDS('Q1_normalized_ensembl.rds')

group_by = 'seurat_clusters'

avg_log2FC_gte = 0.25

p_val_adj_lte = 0.05

min_pct = 0.1

species = 'human'

gene_id = 'ensembl'

lr_network = read.csv('../scSignalMap/inst/extdata/lr_network.csv', header=T)
lr_pairs = na.omit(lr_network[,c(paste('ligand',species,gene_id,sep='_'),paste('receptor',species,gene_id,sep='_'))])


