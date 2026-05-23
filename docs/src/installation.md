---
layout: default
title: Installation
nav_order: 3
---
# Requirements and dependencies

The scSignalMap package is designed to work alongside [Seurat](https://satijalab.org/seurat), and identifies cell-cell interactions from scRNA-Seq data stored in Seurat objects. It is suggested that the latest version of Seurat be installed. The package has been tested on Seurat version 5.5.0.

Installation of scSignalMap from github is facilitated by the [remotes](https://cran.r-project.org/web/packages/remotes/index.html) package. 
scSignalMap usage requires R packages:
`scSignalMap, Seurat, enrichR, dplyr, tidyr, stringr, EnsDb.Hsapiens.v86, AnnotationDbi, fastmatch, and data.table` to run the full pipline. Installation can be done using:

```r
install.packages(c('scSignalMap', 'Seurat', 'enrichR', 'dplyr', 'tidyr', 'stringr', 'EnsDb.Hsapiens.v86', 'AnnotationDbi', 'fastmatch', 'data.table'))
```


#  Installation

Installation of scSignalMap in R is accomplished by:

```r
remotes::install_github("plaisier-lab/scSignalMap/scSignalMap")
```

## Dockerfile and image

For those who are interested we also provide a Dockerfile and image that include all dependencies.
Command to pull the image down:

```r
docker pull cplaisier/quadculture
```
Command to run the docker image. Note that the should be replaced with the path to your files that you want to be mounted onto the docker instance. The files can then be found in /files on the instance and locally on your computer in the path specified.

```r
docker run -it -v '<replace with the location for your files>:/files' cplaisier/quadculture
```

