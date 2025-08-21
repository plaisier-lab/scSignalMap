#' Capturing signaling pathways using scRNA-seq data
#'
#' [Write a description]
#'
#'
#' @param seurat_obj: a seurat object must be supplied to classify, no default
#' @param group_by: a meta_data column that will be used to split the cells, no default
#' @param avg_log2FC_gte: 
#' @param p_val_adj_lte: 
#' @param species: from which species did the samples originate, either 'human' or 'mouse', defaults to 'human'
#' @param gene_id: what type of gene ID is used, either 'ensembl' or 'symbol', defaults to 'ensembl'
#' @return data.frame with putative interactions
#' @export
MapInteractions = function(seurat_obj, group_by, avg_log2FC_gte = 0.25, p_val_adj_lte = 0.05, min_pct = 0.1, species='human', gene_id='ensembl') {
    cat('Running scSignalMap:\n')

    ## Step 1. Load up the data

    cat('  Loading data...\n')
    # Load MultiNicheNet ligand receptor interactions
    lr_network = read.csv(system.file('extdata', 'lr_network.csv', package='scSignalMap'), header=TRUE)
    lr_pairs = na.omit(lr_network[,c(paste('ligand',species,gene_id,sep='_'),paste('receptor',species,gene_id,sep='_'))])
    if(gene_id!='symbol') {
        lig_convert = lr_network[,c(paste('ligand',species,'symbol',sep='_'))]
        names(lig_convert) = lr_network[,c(paste('ligand',species,gene_id,sep='_'))]
        rec_convert = lr_network[,c(paste('receptor',species,'symbol',sep='_'))]
        names(rec_convert) = lr_network[,c(paste('receptor',species,gene_id,sep='_'))]
    }

    # Load up secreted ligands
    secreted = read.csv(system.file('extdata', 'secreted.csv', package='scSignalMap'), header=TRUE)
    secreted_ligands = na.omit(secreted[,paste(species,gene_id,sep='_')])


    ## Step 2. Identify marker genes
    cat('  Identifying marker genes...\n')
    
    putative_markers = FindAllMarkers(seurat_obj, group.by=group_by, min.pct=min_pct, verbose=T)
    markers = list()
    for(cluster1 in sort(unique(seurat_obj@meta.data[,group_by]))) {
        markers[[cluster1]] = putative_markers %>% filter(cluster==cluster1 & avg_log2FC>=avg_log2FC_gte & p_val_adj<=p_val_adj_lte) %>% pull(name='gene')
        #cat(paste0('    ',cluster1,': ',length(markers[[cluster1]]),'\n'))
    }


    ## Step 3. Precompute
    cat('  Precomputing...\n')
    
    # Get a vector of marker genes   
    allgenes = rownames(seurat_obj@assays$RNA)[which(rowSums(as.matrix(seurat_obj@assays$RNA$counts))>0)]

    # Subset ligand_receptor list to those in allgenes
    lr_pairs = lr_pairs[lr_pairs[,1] %fin% allgenes,]
    lr_pairs = lr_pairs[lr_pairs[,2] %fin% allgenes,]

    # Split up the seurat object by groups
    seurat_obj_split = SplitObject(seurat_obj, split.by=group_by)

    # Build a data.table with statistics
    all_dt = rbindlist(lapply(sort(unique(seurat_obj@meta.data[, group_by])), function(clust1) {
                 mat = as.matrix(seurat_obj_split[[clust1]]@assays$RNA@layers$counts)
                 genes = rownames(seurat_obj_split[[clust1]]@assays$RNA)
                 n_cells = ncol(mat)

                 data.table(
                     clust1    = clust1,
                     gene      = genes,
                     counts    = rowSums(mat),
                     perc_gt_0  = rowSums(mat > 0)  / n_cells,
                     perc_gte_3 = rowSums(mat >= 3) / n_cells,
                     perc_gte_10= rowSums(mat >= 10)/ n_cells,
                     avg_exp   = rowMeans(mat)
                 )
             }), use.names = TRUE)


    ## Step 4. Integrate ligand receptor pair with expression data
    
    cat('  L-R and S-R pairs...\n')
    
    
    # All sender/receiver combinations
    clusts = sort(unique(seurat_obj@meta.data[[group_by]]))
    clust_dt = CJ(Sender = clusts, Receiver = clusts)  # Cartesian product

    # Ligand–Receptor pairs table
    if (gene_id == "symbol") {
      lr_dt = data.table(
        Ligand   = lr_pairs[, 1],
        Receptor = lr_pairs[, 2]
      )
    } else {
      lr_dt = data.table(
        Ligand          = lr_pairs[, 1],
        Ligand_Symbol   = lig_convert[lr_pairs[, 1]],
        Receptor        = lr_pairs[, 2],
        Receptor_Symbol = rec_convert[lr_pairs[, 2]]
      )
    }

    # Cartesian product of LR pairs × cluster pairs
    pairs_data = lr_dt[, cbind(.SD, clust_dt), by = seq_len(nrow(lr_dt))][, seq_len := NULL]

    cat(paste0('    Rows = ',nrow(pairs_data),'\n'))

    cat('  Integrating data...\n')
    steps = 13
    pb = txtProgressBar(min = 0, max = steps, style = 3)
    i = 0
    pairs_data[,'Ligand_Counts' := all_dt[pairs_data, on = .(clust1=Sender, gene=Ligand), 'counts']]
    i = i + 1
    setTxtProgressBar(pb, i)
    pairs_data[,'Ligand_gte_3' := all_dt[pairs_data, on = .(clust1=Sender, gene=Ligand), 'perc_gte_3']]
    i = i + 1
    setTxtProgressBar(pb, i)
    pairs_data[,'Ligand_gte_10' := all_dt[pairs_data, on = .(clust1=Sender, gene=Ligand), 'perc_gte_10']]
    i = i + 1
    setTxtProgressBar(pb, i)
    pairs_data[,'Ligand_Cells_Exp' := all_dt[pairs_data, on = .(clust1=Sender, gene=Ligand), 'perc_gt_0']]
    i = i + 1
    setTxtProgressBar(pb, i)
    pairs_data[,'Ligand_Avg_Exp' := all_dt[pairs_data, on = .(clust1=Sender, gene=Ligand), 'avg_exp']]
    i = i + 1
    setTxtProgressBar(pb, i)
    pairs_data[,'Ligand_Cluster_Marker' := mapply(function(sender, ligand) { ligand %in% markers[[sender]] }, Sender, Ligand)]
    i = i + 1
    setTxtProgressBar(pb, i)
    pairs_data[,'Ligand_secreted'] = sapply(pairs_data[,'Ligand'], function(x) { x %fin% secreted_ligands })
    i = i + 1
    setTxtProgressBar(pb, i)
    pairs_data[,'Receptor_Counts' := all_dt[pairs_data, on = .(clust1=Receiver, gene=Receptor), 'counts']]
    i = i + 1
    setTxtProgressBar(pb, i)
    pairs_data[,'Receptor_gte_3' := all_dt[pairs_data, on = .(clust1=Receiver, gene=Receptor), 'perc_gte_3']]
    i = i + 1
    setTxtProgressBar(pb, i)
    pairs_data[,'Receptor_gte_10' := all_dt[pairs_data, on = .(clust1=Receiver, gene=Receptor), 'perc_gte_10']]
    i = i + 1
    setTxtProgressBar(pb, i)
    pairs_data[,'Receptor_Cells_Exp' := all_dt[pairs_data, on = .(clust1=Receiver, gene=Receptor), 'perc_gt_0']]
    i = i + 1
    setTxtProgressBar(pb, i)
    pairs_data[,'Receptor_Avg_Exp' := all_dt[pairs_data, on = .(clust1=Receiver, gene=Receptor), 'avg_exp']]
    i = i + 1
    setTxtProgressBar(pb, i)
    pairs_data[,'Receptor_Cluster_Marker' := mapply(function(receiver, receptor) { receptor %in% markers[[receiver]] }, Receiver, Receptor)]
    i = i + 1
    setTxtProgressBar(pb, i)
    close(pb)
    cat('Done.\n')

    return(pairs_data)
}

