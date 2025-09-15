# scSignalMap: scRNA-Seq Ligand–Receptor Signaling Analysis in R
 > scSignalmap performs a **full ligand–receptor signaling analysis** for single-cell RNA sequencing data. Utilizing user-defined parameters, it performs a streamlined workflow including: 

 - Ligand–receptor interaction mapping
 - Differential expression analysis
 - Receptor filtering
 - Receptor–ligand intersection 
 - Pathway enrichment analysis

 ## Table of Contents

 ## Requirements and Dependencies
scSignalMap usage requires R packages:
`scSignalMap, Seurat, enrichR, dplyr, tidyr, stringr, org.Hs.eg.db, AnnotationDbi, fastmatch, and data.table` to run the full pipline. 
 ## Installation
Installation of scSignalMap in R is accomplished by:

`remotes::install_github("plaisier-lab/scSignalMap/scSignalMap")`

## Dockerfile and Image

## Running scSignalMap
>The scSignalMap pipeline is performed through one function, `run_full_scSignalMap_pipeline()` and an optional function, `create_master_list`, that creates a master list of scSignalMap outputs.

Below is an example of inputs for `run_full_scSignalMap_pipeline`:

```
results = run_full_scSignalMap_pipeline(
workingdir = "/files", 
seurat_obj = 'seurat_objects/MN_big_int.rds', 
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
enrichr_databases = c("BioCarta_2016", 
		    "GO_Biological_Process_2025", 
		    "KEGG_2021_Human", 	
		    "NCI-Nature_2016", 	
		    "WikiPathways_2024_Human"), 
adj_p_val_method = "BH")
```
Below is an example of inputs for `create_master_list`:



## Accessing Outputs
Accessing results from `run_full_scSignalMap_pipeline` can be accomplished using:
`write.csv()`
## Maintainers
