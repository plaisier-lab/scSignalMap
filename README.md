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
Seurat, dplyr, enrichR, org.Hs.eg.db, data.table, 
 ## Installation
Installation of scSignalMap in R is accomplished by:

`remotes::install_github("plaisier-lab/scSignalMap/scSignalMap")`

## Dockerfile and Image

## Running scSignalMap
The scSignalMap pipeline is performed through one function, `run_full_scSignalMap_pipeline()` and an optional function that creates a master list of scSignalMap outputs.

## Accessing Outputs

## Maintainers
