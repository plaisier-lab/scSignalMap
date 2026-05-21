---
layout: default
title: Running scSignalMap
nav_order: 4
---
# Running scSignalMap
>The scSignalMap pipeline is performed through one function, `run_full_scSignalMap_pipeline()` and an optional function, `create_master_interaction_list`, that creates a master list of scSignalMap outputs.

Below is an example of inputs for the `run_full_scSignalMap_pipeline` function:

```r
    results = run_full_scSignalMap_pipeline(
        seurat_obj = seurat1, 
        prep_SCT = TRUE, 
        cond_column = 'sample', 
        cond_name1 = paste0(GB3_BBB'), 
        cond_name2 = paste0('GB3_MVN'), 
        celltype_column = 'celltype', 
        celltype_name = 'Tumor', 
        sender_celltypes = c('Astrocytes', 'Endothelial_Cells', 'Pericytes'), 
        receiver_celltypes = 'Tumor', 
        secreted_lig = TRUE, 
        FC_cutoff = 0.3, 
        adj_p_val_cutoff = 0.05, 
        enrichr_databases = c('BioCarta_2016', 'GO_Biological_Process_2025', 'KEGG_2021_Human', 'NCI-Nature_2016', 'WikiPathways_2024_Human'), 
        adj_p_val_method = 'BH')
```
