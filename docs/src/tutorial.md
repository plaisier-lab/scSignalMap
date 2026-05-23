---
layout: default
title: Tutorial Using Glioblastoma (GBM) scRNA-Seq Datasets
nav_order: 6
---

# Input Information
The scSignalMap requires a thoroughly quality-controlled Seurat
object. 

An example QC pipeline used on the GBM datasets is available
[here].

## GBM Input Data

The GBM dataset used to test the scSignalMap is
available
[here].
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

