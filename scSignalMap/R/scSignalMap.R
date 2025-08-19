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
    
    putative_markers = FindAllMarkers(seurat_obj, group.by=group_by, min.pct=min_pct, verbose=FALSE)
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

    # Prepare lists to hold precomputed data
    counts = list()
    perc_gt0 = list()
    perc_gte3 = list()
    perc_gte10 = list()
    avg_exp = list()

    # Split up the seurat object by groups
    seurat_obj_split = SplitObject(seurat_obj, split.by=group_by)

    # Iterate through each cluster
    for(clust1 in sort(unique(seurat_obj@meta.data[,group_by]))) {
        #cat(paste0('    ',clust1,'\n'))
        # Precompute counts per cluster
        cnts = rowSums(as.matrix(seurat_obj_split[[clust1]]@assays$RNA@layers$counts))
        names(cnts) = rownames(seurat_obj_split[[clust1]]@assays$RNA)
        counts[[clust1]] = cnts
        
        # Precompute percentage of cells with at least 1 transcript per cluster
        rowSums_gt0 = rowSums(as.matrix(seurat_obj_split[[clust1]]@assays$RNA@layers$counts)>0)/ncol(seurat_obj_split[[clust1]])
        names(rowSums_gt0) = rownames(seurat_obj_split[[clust1]]@assays$RNA)
        perc_gt0[[clust1]] = rowSums_gt0

        # Precompute number of cells with at least 3 transcript per cluster
        rowSums_gte3 = rowSums(as.matrix(seurat_obj_split[[clust1]]@assays$RNA@layers$counts)>=3)/ncol(seurat_obj_split[[clust1]])
        names(rowSums_gte3) = rownames(seurat_obj_split[[clust1]]@assays$RNA)
        perc_gte3[[clust1]] = rowSums_gte3
        
        # Precompute number of cells with at least 3 transcript per cluster
        rowSums_gte10 = rowSums(as.matrix(seurat_obj_split[[clust1]]@assays$RNA@layers$counts)>=10)/ncol(seurat_obj_split[[clust1]])
        names(rowSums_gte10) = rownames(seurat_obj_split[[clust1]]@assays$RNA)
        perc_gte10[[clust1]] = rowSums_gte10
        
        # Precompute average expression per cluster
        avg_clust1 = rowMeans(as.matrix(seurat_obj_split[[clust1]]@assays$RNA@layers$counts))
        names(avg_clust1) = rownames(seurat_obj_split[[clust1]]@assays$RNA)
        avg_exp[[clust1]] = avg_clust1
    }


    ## Step 4. Integrate ligand receptor pair with expression data
    
    cat('  Integration...\n')
    
    
    # Iterate through ligand receptor pairs
    t0 = proc.time()[3]
    i = 1
    j = 1
    pb = txtProgressBar(min = 0, max = nrow(lr_pairs), style=3, width=50, char= '=')
    pairs_data = list()
    for(pair1 in 1:nrow(lr_pairs)) {
        lig1 = lr_pairs[pair1,1]
        rec1 = lr_pairs[pair1,2]
        # Make sure ligand receptor pair is in expression data
        #cat(paste0('    ',lig1,'->',rec1,'\n'))
        
        # Iterate through sender cell types
        for(clust1 in sort(unique(seurat_obj@meta.data[,group_by]))) {
            if(gene_id=='symbol') {
                tmp1 = c(lig1, rec1, clust1)
            } else {
                tmp1 = c(lig1, lig_convert[lig1], rec1, rec_convert[rec1], clust1)
            }
            tmp2 = c(tmp1, counts[[clust1]][lig1], perc_gte3[[clust1]][lig1], perc_gte10[[clust1]][lig1], perc_gt0[[clust1]][lig1], avg_exp[[clust1]][lig1])

            # Iterate through reciever cell types
            for(clust2 in sort(unique(seurat_obj@meta.data[,group_by]))) {
                # Row bind the data into the matrix
                tmp3 = c(tmp2, clust2, counts[[clust2]][rec1], perc_gte3[[clust2]][rec1], perc_gte10[[clust2]][rec1], perc_gt0[[clust2]][rec1], avg_exp[[clust2]][rec1])
                pairs_data[[i]] = tmp3
                i = i + 1
            }
        }
        setTxtProgressBar(pb, j)
        j = j + 1
    }
    
    # Prepare a matrix to hold the data
    pairs_data = data.frame(do.call(rbind, pairs_data))
    if(gene_id=='symbol') {
        colnames(pairs_data) = c('Ligand','Receptor','Sender','Receiver')
    } else {
        colnames(pairs_data) = c('Ligand','Ligand_Symbol','Receptor','Receptor_Symbol','Sender','Receiver')
    }
    t0_1 = proc.time()[3]
    print(paste0('Total time: ',t0_1-t0))


    t1 = proc.time()[3]
    pairs_data[,'Ligand_Cluster_Marker'] = sapply(1:nrow(pairs_data), function(x) { pairs_data[x,'Ligand'] %fin% markers[[pairs_data[x,'Sender']]]})
    t2 = proc.time()[3]
    print(paste0('Ligand_Cluster_Marker: ',t2-t1))
    t1 = proc.time()[3]
    pairs_data[,'Ligand_secreted'] = pairs_data[,'Ligand'] %fin% secreted_ligands
t2 = proc.time()[3]
    print(paste0('Ligand_secreted: ',t2-t1))
    t1 = proc.time()[3]
    pairs_data[,'Receptor_Cluster_Marker'] = sapply(1:nrow(pairs_data), function(x) { pairs_data[x,'Receptor'] %fin% markers[[pairs_data[x,'Receiver']]]})
t2 = proc.time()[3]
    print(paste0('Receptor_Cluster_Marker: ',t2-t1))
    cat('Done.\n')

    return(pairs_data)
}


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
MapInteractions_vec = function(seurat_obj, group_by, avg_log2FC_gte = 0.25, p_val_adj_lte = 0.05, min_pct = 0.1, species='human', gene_id='ensembl') {
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
    
    putative_markers = FindAllMarkers(seurat_obj, group.by=group_by, min.pct=min_pct, verbose=FALSE)
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
    pairs_data[,'Ligand_Counts' := all_dt[pairs_data, on = .(clust1=Sender, gene=Ligand), 'counts']]
    pairs_data[,'Lig_gte_3' := all_dt[pairs_data, on = .(clust1=Sender, gene=Ligand), 'perc_gte_3']]
    pairs_data[,'Lig_gte_10' := all_dt[pairs_data, on = .(clust1=Sender, gene=Ligand), 'perc_gte_10']]
    pairs_data[,'Ligand_Cells_Exp' := all_dt[pairs_data, on = .(clust1=Sender, gene=Ligand), 'perc_gt_0']]
    pairs_data[,'Ligand_Avg_Exp' := all_dt[pairs_data, on = .(clust1=Sender, gene=Ligand), 'avg_exp']]
    pairs_data[,'Ligand_Cluster_Markter' := mapply(function(sender, ligand) { ligand %in% markers[[sender]] }, Sender, Ligand)]
    pairs_data[,'Ligand_secreted'] = pairs_data[,'Ligand'] %fin% secreted_ligands
    pairs_data[,'Receptor_Counts' := all_dt[pairs_data, on = .(clust1=Receiver, gene=Receptor), 'counts']]
    pairs_data[,'Rec_gte_3' := all_dt[pairs_data, on = .(clust1=Receiver, gene=Receptor), 'perc_gte_3']]
    pairs_data[,'Rec_gte_10' := all_dt[pairs_data, on = .(clust1=Receiver, gene=Receptor), 'perc_gte_10']]
    pairs_data[,'Receptor_Cells_Exp' := all_dt[pairs_data, on = .(clust1=Receiver, gene=Receptor), 'perc_gt_0']]
    pairs_data[,'Receptor_Avg_Exp' := all_dt[pairs_data, on = .(clust1=Receiver, gene=Receptor), 'avg_exp']]
    pairs_data[,'Receptor_Cluster_Marker' := mapply(function(receiver, receptor) { receptor %in% markers[[receiver]] }, Receiver, Receptor)]
    cat('Done.\n')

    return(pairs_data)
}

