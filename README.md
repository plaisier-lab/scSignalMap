# scSignalMap: scRNA-Seq Ligand–Receptor Signaling Analysis in R
 > scSignalmap performs a **full ligand–receptor signaling analysis** for single-cell RNA sequencing data. Utilizing user-defined parameters, it performs a streamlined workflow including: 

 - Ligand–receptor interaction mapping
 - Differential expression analysis
 - Receptor filtering
 - Receptor–ligand intersection 
 - Pathway enrichment analysis

 ## Table of Contents
- Requirements and Dependencies
- Installation
- Dockerfile and Image
- Tutorials
  - Running scSignalMap
  - Accessing Outputs
- Maintainers

 ## Requirements and Dependencies
scSignalMap usage requires R packages:
`scSignalMap, Seurat, enrichR, dplyr, tidyr, stringr, org.Hs.eg.db, AnnotationDbi, fastmatch, and data.table` to run the full pipline. 
 ## Installation
Installation of scSignalMap in R is accomplished by:
```r
remotes::install_github("plaisier-lab/scSignalMap/scSignalMap")
```
## Dockerfile and Image

## Running scSignalMap
>The scSignalMap pipeline is performed through one function, `run_full_scSignalMap_pipeline()` and an optional function, `create_master_interaction_list`, that creates a master list of scSignalMap outputs.

Below is an example of inputs for `run_full_scSignalMap_pipeline`:

```r
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
Below is an example of inputs for `create_master_interaction_list`:

```r
master_interaction_list = create_master_interaction_list(
  enrichr_results = results$enrichr_results,
  de_receptors = results$upreg_receptors_filtered_and_compared,
  scSignalMap_data_filtered = results$interactions_filtered)
```

## Accessing Outputs
Accessing results from `run_full_scSignalMap_pipeline` can be accomplished using the following code:
```r
write.csv(results$enrichr_results, "<file_name.csv>")
```
Accessing results from `create_master_interaction_list` can be accomplished using the following code:
```r
write.csv(master_interaction_list, "<file_name.csv>")
```
## Maintainers
For issues or comments, please contact: [Chris Plaiser](mailto:plaisier@asu.edu)

For other great packages from the Plaisier Lab, please check here: [@plaisier_lab](https://github.com/plaisier_lab)
