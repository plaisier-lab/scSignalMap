---
layout: default
title: Tutorial Using Glioblastoma (GBM) scRNA-Seq Datasets
nav_order: 6
---

# Input Information
The scSignalMap requires a thoroughly quality-controlled Seurat
object. 

An example QC pipeline used on the GBM datasets is available
[here](https://github.com/plaisier-lab/scSignalMap/blob/main/test/GBM_QC.R).

## GBM Input Data

The GBM dataset used to test the scSignalMap is
available [here].
This data has been QC'd and normalized using SCTransform following our
best practices described above.

### Load Packages 
```r
library(scSignalMap)
library(Seurat)
library(enrichR)
library(dplyr)
library(tidyr)
library(stringr)
library(EnsDb.Hsapiens.v86)
library(AnnotationDbi)
library(fastmatch)
library(data.table)
```

### Load in Sample Data
```r
seurat_obj = readRDS("GB3_Integrated_Dataset.rds")
```

### Run scSignalMap on data
```r
results = run_full_scSignalMap_pipeline(
	seurat_obj = seurat_obj, 
	prep_SCT = TRUE, 
	cond_column = "orig.ident", 
	cond_name1 = "MN2", 
	cond_name2 = "MN1", 
	celltype_column = "celltype", 
	celltype_name = "GB3_Tumor", 
	sender_celltypes = c("Astrocyte", "HUVEC"), 
	receiver_celltypes = "GB3_Tumor", 
	secreted_lig = TRUE, 
	FC_cutoff = 0.3, 
	adj_p_val_cutoff = 0.05, 
	enrichr_databases = c("BioCarta_2016", "GO_Biological_Process_2025", "KEGG_2021_Human", "NCI-Nature_2016", "WikiPathways_2024_Human"), 
	adj_p_val_method = "BH")
```

### Access Results
```r
#Master List Call:
master_interaction_list = create_master_interaction_list(
  enrichr_results = results$enrichr_results,
  de_receptors = results$de_receptors_filtered_and_compared,
  scSignalMap_data_filtered = results$interactions_filtered
)
write.csv(master_interaction_list, "master_interaction_list.csv", row.names = FALSE)
```

The output from this query can be found [here]
