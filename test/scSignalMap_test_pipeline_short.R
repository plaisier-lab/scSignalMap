#----docker 
#docker run -it -v "/home/jswoodl2/MN_data:/files" --rm cplaisier/quadculture

remotes::install_github('plaisier-lab/scSignalMap/scSignalMap@v1.1')
install.packages("enrichR")
install.packages("fastmatch")
library(scSignalMap)
library(enrichR)
library(Seurat)
library(dplyr)
library(tidyr)
library(stringr)
library(EnsDb.Hsapiens.v86)
library(AnnotationDbi)
library(fastmatch)
library(data.table)


###################################################################
############# ---------- Test Out Analysis ---------- #############
###################################################################

setwd('/files')
seurat1 = readRDS("seurat_objects/MN_big_int.rds")


## Run function and save as .csv (add gene_id='symbol' if needed)
LR_interactions = MapInteractions(seurat1,
                                  'celltype')
write.csv(LR_interactions, "MN_big_int_scSignalMap_interactions.csv")


## Find DE genes between MN2 and MN1 for GB3_Tumor cells
de_cond_celltype = find_markers_btwn_cond_for_celltype(
    seurat_obj = seurat1,
    prep_SCT = TRUE,
    cond_column = "orig.ident",
    cond_name1 = "MN2",
    cond_name2 = "MN1",
    celltype_column = "celltype",
    celltype_name = "GB3_Tumor",
    FC_cutoff = 0.3,
    adj_p_val_cutoff = 0.05)
write.csv(de_cond_celltype, "de_btwn_MN2_vs_MN1_GB3_Tumor.csv", row.names = FALSE)


## Find upregulated receptors
upreg_receptors = find_upreg_receptors(
    de_condition_filtered = de_cond_celltype,
    FC_cutoff = 0.3)
write.csv(upreg_receptors, "upregulated_receptors_btwn_MN2_vs_MN1_GB3_Tumor.csv", row.names = FALSE)

## Filter LR Interactions
interactions_filtered = filter_lr_interactions(
    interactions = LR_interactions,
    sender_celltypes = c("Astrocyte", "HUVEC"),
    receiver_celltypes = "GB3_Tumor",
    secreted_lig = TRUE)


## Intersect upregulated receptors with ligand-receptor pairs
upreg_receptors_filtered_and_compared = intersect_upreg_receptors_with_lr_interactions(
    upreg_receptors = upreg_receptors,
    interactions = interactions_filtered)
write.csv(upreg_receptors_filtered_and_compared, "upreg_receptors_filtered_and_compared.csv", row.names = FALSE)

## Find pathways that are enriched with DE genes that include the upregulated receptor
enrichr_results = find_enriched_pathways(
    seurat_obj = seurat1,
    de_condition_filtered = de_cond_celltype,
    enrichr_databases = c("BioCarta_2016", "GO_Biological_Process_2025", "KEGG_2021_Human", "NCI-Nature_2016", "WikiPathways_2024_Human"),
    adj_p_val_method = "BH",
    adj_p_val_cutoff = 0.05)
write.csv(enrichr_results, "enrichr_results.csv", row.names = FALSE)


## Master List Call 
master_interaction_list = create_master_interaction_list(
  enrichr_results = enrichr_results,
  de_receptors = upreg_receptors_filtered_and_compared,
  scSignalMap_data_filtered = interactions_filtered)
write.csv(master_interaction_list, "12.02.25_master_interaction_list.csv", row.names = FALSE)

