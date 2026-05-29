#--------------------------------
# Set up section / load packages
#--------------------------------

library(dplyr)
library(Seurat)
library(SeuratDisk)
library(keras)
library(patchwork)
library(ggplot2)
library(grid)
library(gridExtra)
library(writexl)
library(data.table)
library(ccAFv2)
library("org.Hs.eg.db")
library(aricode)
library(copykat)

# Set working directory
setwd('/files')

reticulate::use_python('/usr/bin/python3.9')

# Mitochondrial genes as ensembl IDs
mito_genes = read.csv("mito_genes.csv")[['mito']]

# Load ccSeurat phase gene sets
s.genes = cc.genes$s.genes
g2m.genes = cc.genes$g2m.genes
# convert to ensembl IDs
ensembl_s_genes = mapIds(org.Hs.eg.db, keys = s.genes, keytype = "SYMBOL", column="ENSEMBL", multiVals='first')
ensembl_s_genes_2 = na.omit(data.frame(ensembl_s_genes))$ensembl_s_genes
ensembl_g2m_genes = mapIds(org.Hs.eg.db, keys = g2m.genes, keytype = "SYMBOL", column="ENSEMBL", multiVals='first')
ensembl_g2m_genes_2 = na.omit(data.frame(ensembl_g2m_genes))$ensembl_g2m_genes


#-------------
# Functions
#-------------

scQC = function(tag, mito_genes, v1 = 2000, v2 = 150000, h1 = 0.001, h2 = 0.1, data_dir = 'outs', save_dir = 'analysis_output', obj_dir = 'seurat_objects', mt = 'MT-', symbol = F, resolution=0.5) {
  cat('\n',tag,'\n')

  # Create folders
  dir.create(save_dir, showWarnings = FALSE)
  resdir2 = save_dir
  dir.create(obj_dir, showWarnings = FALSE)
  resdir3 = obj_dir
  dir.create('copykat', showWarnings = FALSE)

  #---------------------
  # Load in data
  #---------------------
  gene_column = 1
  gene_id = 'ensembl'
  if(symbol){
    gene_column = 2
    gene_id = 'gene_symbols'
  }
  data = Read10X(data_dir, gene.column=gene_column) 
  cat('Raw data', dim(data)[2], 'cells', dim(data)[1], 'genes \n')
  rownames(data) = gsub("_", "-", rownames(data))

  # Create seurat object
  seurat1 = CreateSeuratObject(counts = data, min.cells = 3, min.features = 200, project=tag)
  cat('Basic filter', dim(seurat1)[2], 'cells', dim(seurat1)[1], 'genes \n')
  if(symbol){
    mito_genes = grep(mt, rownames(seurat1))
  }
  seurat1[['percent.mito']] = PercentageFeatureSet(seurat1, features = mito_genes)/100

  #---------------------
  # Quality control
  #---------------------
  cat('Quality control \n')

  # Quality control plots for choosing cutoffs
  pdf(file.path(resdir2, paste0(tag, '_QC_plot_to_choose_cutoffs.pdf')))
  plot(seurat1@meta.data$nCount_RNA, seurat1@meta.data$percent.mito,
       xlab = 'nCount_RNA', ylab = 'percent.mito', pch = 20)
  abline(v = v1, col = 'red', lwd =3, lty =2)
  text(v1,0,as.character(v1), cex = 0.75, pos = 1)
  abline(v = v2, col = 'red', lwd =3, lty =2)
  text(v2,0,as.character(v2), cex = 0.75, pos = 1)
  abline(h = h1 , col = 'red', lwd =3, lty =2)
  text(as.character(v2+10000),h1,as.character(h1), cex = 0.75, pos = 3)
  abline (h = h2, col = 'red', lwd =3, lty =2)
  text(as.character(v2+10000),h2,as.character(h2), cex = 0.75, pos = 3)
  print(VlnPlot(seurat1, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3))
  dev.off()

  # Quality control filtering
  keep.detect = which(seurat1@meta.data$percent.mito < h2 & seurat1@meta.data$percent.mito > h1 & seurat1@meta.data$nCount_RNA < v2 & seurat1@meta.data$nCount_RNA > v1)
  seurat1 = subset(seurat1, cells=colnames(seurat1)[keep.detect])


  # Run copykat
   if (!file.exists(file.path('copykat',paste0(tag, '_copykat_prediction.txt')))) {
     exp.rawdata = as.matrix(seurat1@assays$RNA@layers$counts)
     dimnames(exp.rawdata) = dimnames(seurat1@assays$RNA)
     copykat_res = copykat(rawmat=exp.rawdata, id.type="E", ngene.chr=5, win.size=25, KS.cut=0.1, sam.name=paste0('copykat/',tag), distance="euclidean", norm.cell.names="", output.seg="FALSE", plot.genes="TRUE", genome="hg20",n.cores=12)
   }
  copykat_pred = read.csv(file.path('copykat', paste0(tag, '_copykat_prediction.txt')), header = TRUE, sep = '\t')

  copykat_pred = copykat_pred[match(colnames(seurat1), copykat_pred$cell.names), ]
  seurat1$copykat = copykat_pred$copykat.pred

   # Save out filtered object
   cat('Filtered to', dim(seurat1)[2], 'cells', dim(seurat1)[1], 'genes \n')
   saveRDS(seurat1, file.path(resdir3, paste0(tag, '_filtered_', paste0(gene_id), '.rds')))

  seurat2 = seurat1

  #---------------------------------
  # Normalization with sctransform
  #---------------------------------
  cat('Normalization \n')
  seurat2 = SCTransform(seurat2, verbose = FALSE)
  cat('Normalized genes:', dim(seurat2@assays$SCT@data)[1], 'features,', length(seurat2@assays$SCT@var.features), 'highly variable genes \n')

  # Classify with ccSeurat and save out as csv
  if(symbol){
    seurat2 = CellCycleScoring(object=seurat2, s.features=s.genes, g2m.features=g2m.genes, set.ident=FALSE)
    write.csv(seurat2$Phase, file.path(resdir2, paste0(tag, '_ccSeurat_calls.csv')))
  }
  seurat2 = PredictCellCycle(seurat2, include_g0=T)

  # UMAP
  seurat2 = RunPCA(seurat2, dims = 1:30, verbose=FALSE)
  seurat2 = FindNeighbors(seurat2, dims = 1:30, verbose=FALSE)
  seurat2 = FindClusters(seurat2, verbose=FALSE, resolution = resolution)
  seurat2 = RunUMAP(seurat2, dims=1:30, verbose=FALSE)
  pdf(file.path(resdir2, paste0(tag, '_UMAPs.pdf')))
  d1 = DimPlot(seurat2, reduction = "umap", label=F, group.by="seurat_clusters", raster = FALSE) + ggtitle("seurat_clusters")
  print(d1)
  d2 = DimPlot(seurat2, reduction = "umap", label=F, group.by='copykat', raster = FALSE) + ggtitle("seurat_clusters")
  print(d2)

  #Use gene symbols:
  goi1 = c('PECAM1','S100B','CDKN2A')
  features1 = unlist(mapIds(org.Hs.eg.db, keys = goi1, keytype = "SYMBOL", column="ENSEMBL", multiVals='first'))
  
  for(i in c(1:length(features1))){
    f1 = FeaturePlot(seurat2, reduction = "umap", features=c(features1[i]), raster = FALSE) + ggtitle(names(goi1[i]))
    print(f1)
  }
  f2 = FeaturePlot(seurat2, reduction = "umap", features=c('ENSG00000261371'), raster = FALSE) + ggtitle('PECAM1')
  print(f2)
  f3 = FeaturePlot(seurat2, reduction = "umap", features=c('ENSG00000160307'), raster = FALSE) + ggtitle('S100B')
  print(f3)
  dev.off()

  cluster_markers= FindAllMarkers(seurat2, logfc.threshold = 0.25, only.pos = TRUE)
  cluster_markers_genes = mapIds(org.Hs.eg.db, keys = cluster_markers$gene, keytype = "ENSEMBL", column="SYMBOL", multiVals='first')
  cluster_markers$gene = cluster_markers_genes
  write.csv(cluster_markers, file.path('analysis_output', paste0(tag,"_scTransform_Markers_together.csv")))
  cat('Saving normalized RDS object \n')
  saveRDS(seurat2, file.path(resdir3, paste0(tag, '_normalized_', paste0(gene_id), '.rds')))
  return(seurat2)
}


#--------------------------------------
# Quality control & data preparation
#--------------------------------------

#MN Data
pip_ensembl = list()

# Monoculture 
pip_ensembl[['MN1']] = scQC(data_dir = 'MN1_monoculture/filtered_feature_bc_matrix', tag = 'MN1', mito_genes = mito_genes, v1 = 10000, v2 = 180000, h1 = 0.01, h2 = 0.18, resolution=0.02, save_dir = 'analysis_output', obj_dir = 'seurat_objects', mt = 'MT-', symbol = F)

# Triculture
pip_ensembl[['MN2']] = scQC(data_dir = 'MN2_triculture/filtered_feature_bc_matrix', tag = 'MN2', mito_genes = mito_genes, v1 = 10000, v2 = 190000, h1 = 0.01, h2 = 0.16, resolution=0.02)

#--------------------------------------
# Merge Monoculture and Triculture: 
#--------------------------------------
MN_big = merge(pip_ensembl[['MN1']], y = c(pip_ensembl[['MN2']]), add.cell.ids=c('MN1','MN2'))

MN_big = SCTransform(MN_big, verbose = FALSE, return.only.var.genes = FALSE)
MN_big = RunPCA(MN_big, dims = 1:30, verbose=FALSE)
MN_big_int = IntegrateLayers(object = MN_big, method = CCAIntegration, orig.reduction = "pca", new.reduction = "integrated.cca", normalization.method='SCT', verbose = FALSE)
MN_big_int = JoinLayers(MN_big_int, assay='RNA')

MN_big_int = RunPCA(MN_big_int, dims = 1:30, verbose=FALSE)
MN_big_int = FindNeighbors(MN_big_int, dims = 1:30, verbose=FALSE)
MN_big_int = FindClusters(MN_big_int, verbose=FALSE, resolution = 0.02)
MN_big_int = RunUMAP(MN_big_int, dims=1:30, verbose=FALSE)

pdf(file.path('analysis_output' ,paste0('MN_big_UMAPs.pdf')))
d1 = DimPlot(MN_big_int, reduction = "umap", label=F, group.by="seurat_clusters", raster = FALSE) + ggtitle("seurat_clusters")
print(d1)
d2 = DimPlot(MN_big_int, reduction = "umap", label=F, group.by='orig.ident', raster = FALSE) + ggtitle("orig.ident")
print(d2)
d3 = DimPlot.ccAFv2(MN_big_int)
print(d3)
f1 = FeaturePlot(MN_big_int, reduction = "umap", features=c('ENSG00000261371'), raster = FALSE) + ggtitle('PECAM1')
print(f1)
f3 = FeaturePlot(MN_big_int, reduction = "umap", features=c('ENSG00000147889'), raster = FALSE) + ggtitle('CDKN2A')
print(f3)
f4 = FeaturePlot(MN_big_int, reduction = "umap", features=c('ENSG00000160307'), raster = FALSE) + ggtitle('S100B')
print(f4)
dev.off()

MN_big_int = PrepSCTFindMarkers(MN_big_int)
cluster_markers= FindAllMarkers(MN_big_int, logfc.threshold = 0.25, only.pos = TRUE) 
cluster_markers_genes = mapIds(org.Hs.eg.db, keys = cluster_markers$gene, keytype = "ENSEMBL", column="SYMBOL", multiVals='first')
cluster_markers$gene = cluster_markers_genes

write.csv(cluster_markers, file.path('analysis_output', paste0('MN_big_int_scTransform_Markers_together.csv')))

# Rename clusters 
new.cluster.ids = c(
  "0" = "HUVEC",
  "1" = "GB3_Tumor",
  "2" = "Astrocyte")

# Rename identities
MN_big_int = RenameIdents(MN_big_int, new.cluster.ids)

# Store renamed identities as metadata
MN_big_int$celltype = Idents(MN_big_int)

# Verify mapping
table(MN_big_int$celltype)

pdf(file.path('analysis_output' ,paste0('MN_big_int_UMAPs_renamed.pdf')))
p1 = DimPlot(MN_big_int, reduction = "umap", label = F, pt.size = 0.5) + NoLegend()
p1 = LabelClusters(p1, id='ident', fontface='bold', color='black')
print(p1)
DimPlot(MN_big_int, reduction = "umap", label = F, pt.size = 0.5)
dev.off()

#Save integrated object
saveRDS(MN_big_int, file = "seurat_objects/GB3_Integrated_Dataset.rds")



