# docker run -it -v '/home/cplaisier/scSignalMap:/files' cplaisier/quadculture

## Install
# remotes::install_github('plaisier-lab/scSignalMap/scSignalMap')
# install.packages('enrichR')

library(Seurat)
library(dplyr)
library(fastmatch)
library(data.table)
library(scSignalMap)
library(EnsDb.Mmusculus.v79)
library(enrichR)

setwd('/files/test')

seurat_obj = readRDS('all_combined_sct.rds')

result_scSignalMap_pipeline_Endo_VSMCFibro = run_full_scSignalMap_pipeline(seurat_obj = seurat_obj, prep_SCT = TRUE, cond_column = 'source', 
                                                                             cond_name1 = 'MoNP', cond_name2 = 'MoNP_VP', celltype_column = 'cell_types', 
                                                                             celltype_name = 'VSMC/Fibro.', sender_celltypes = 'Endo.', receiver_celltypes = 'VSMC/Fibro.', 
                                                                             secreted_lig = TRUE, FC_cutoff = 0.3, adj_p_val_cutoff = 0.05, enrichr_databases = 
                                                                               c("BioCarta_2016", "GO_Biological_Process_2025", "KEGG_2019_Mouse", "NCI-Nature_2016", "WikiPathways_2024_Mouse"), 
                                                                             adj_p_val_method = "BH", ensdb = 'EnsDb.Mmusculus.v79', species='mouse')


result_scSignalMap_pipeline_Endo_Mo = run_full_scSignalMap_pipeline(seurat_obj = seurat_obj, prep_SCT = TRUE, cond_column = 'source', cond_name1 = 'MoNP', 
                                                                      cond_name2 = 'MoNP_VP', celltype_column = 'cell_types', celltype_name = 'Mo', sender_celltypes = 'Endo.', 
                                                                      receiver_celltypes = 'Mo', secreted_lig = TRUE, FC_cutoff = 0.3, adj_p_val_cutoff = 0.05, enrichr_databases = 
                                                                        c("BioCarta_2016", "GO_Biological_Process_2025", "KEGG_2019_Mouse", "NCI-Nature_2016", "WikiPathways_2024_Mouse"), 
                                                                      adj_p_val_method = "BH", ensdb = 'EnsDb.Mmusculus.v79', species='mouse')

