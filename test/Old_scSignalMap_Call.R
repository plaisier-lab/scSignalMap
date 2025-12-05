# docker run -it -v '/home/cplaisier/scSignalMap:/files' cplaisier/ccafv2_s5_extra

## Install
# install.packages('fastmatch')
# remotes::install_github('plaisier-lab/scSignalMap/scSignalMap')

library(Seurat)
library(dplyr)
library(fastmatch)
library(data.table)
library(scSignalMap)

setwd('/files/test')

seurat_obj = readRDS('Q1_normalized_ensembl.rds')

interactions1 = MapInteractions(seurat_obj, 'seurat_clusters')
write.csv(interactions1, 'Q1_norm_ens_scSignalMap_new.csv')

