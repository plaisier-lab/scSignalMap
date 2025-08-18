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
    putative_markers = FindAllMarkers(seurat_obj, group.by=group_by, min.pct=min_pct)
    markers = list()
    for(cluster1 in sort(unique(seurat_obj@meta.data[,group_by]))) {
        markers[[cluster1]] = putative_markers %>% filter(cluster==cluster1 & avg_log2FC>=avg_log2FC_gte & p_val_adj<=p_val_adj_lte) %>% pull(name='gene')
        cat(paste0('    ',cluster1,': ',length(markers[[cluster1]]),'\n'))
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
        cat(paste0('    ',clust1,'\n'))
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
    
    cat('  Integrating data..\n')
    
    # Prepare a matrix to hold the data
    if(gene_id=='symbol') {
        pairs_data = matrix(nrow=0, ncol=17)
        colnames(pairs_data) = c('Ligand','Receptor', 'Sender', 'Ligand_Counts', 'Lig_gte_3', 'Lig_gte_10', 'Ligand_Cells_Exp', 'Ligand_Avg_Exp', 'Ligand_Cluster_Marker', 'Lig_secreted', 'Receiver', 'Receptor_Counts', 'Rec_gte_3', 'Rec_gte_10', 'Receptor_Cells_Exp', 'Receptor_Avg_Exp', 'Receptor_Cluster_Marker')
    } else {
        pairs_data = matrix(nrow=0, ncol=19)
        colnames(pairs_data) = c('Ligand','Ligand Symbol','Receptor','Receptor Symbol', 'Sender', 'Ligand_Counts', 'Lig_gte_3', 'Lig_gte_10', 'Ligand_Cells_Exp', 'Ligand_Avg_Exp', 'Ligand_Cluster_Marker', 'Lig_secreted', 'Receiver', 'Receptor_Counts', 'Rec_gte_3', 'Rec_gte_10', 'Receptor_Cells_Exp', 'Receptor_Avg_Exp', 'Receptor_Cluster_Marker')
    }

    ligMarker = FALSE
    ligSec = FALSE
    recMarker = FALSE
    # Iterate through ligand receptor pairs
    for(pair1 in 1:nrow(lr_pairs)) {
        lig1 = lr_pairs[pair1,1]
        rec1 = lr_pairs[pair1,2]
        # Make sure ligand receptor pair is in expression data
        cat(paste0('    ',lig1,'->',rec1,'\n'))
        
        # Iterate through sender cell types
        for(clust1 in sort(unique(seurat_obj@meta.data[,group_by]))) {
            # Add if ligand is DEG
            #ligMarker = lig1 %fin% markers[[clust1]]
            
            # Add if ligand is secreted
            #ligSec = lig1 %fin% secreted_ligands
            if(gene_id=='symbol') {
                tmp1 = c(lig1, rec1, clust1, counts[[clust1]][lig1], perc_gte3[[clust1]][lig1], perc_gte10[[clust1]][lig1], perc_gt0[[clust1]][lig1], avg_exp[[clust1]][lig1], ligMarker, ligSec)
            } else {
                tmp1 = c(lig1, lig_convert[lig1], rec1, rec_convert[rec1], clust1, counts[[clust1]][lig1], perc_gte3[[clust1]][lig1], perc_gte10[[clust1]][lig1], perc_gt0[[clust1]][lig1], avg_exp[[clust1]][lig1], ligMarker, ligSec)
            }
            
            # Iterate through reciever cell types
            for(clust2 in sort(unique(seurat_obj@meta.data[,group_by]))) {
                # Add if receptor is DEG
                #recMarker = rec1 %fin% markers[[clust2]]

                # Row bind the data into the matrix
                tmp2 = c(tmp1, clust2, counts[[clust2]][rec1], perc_gte3[[clust2]][rec1], perc_gte10[[clust2]][rec1], perc_gt0[[clust2]][rec1], avg_exp[[clust2]][rec1], recMarker)
                pairs_data = rbind(pairs_data, tmp2)
            }
        }
    }

    cat('Done\n')
    return(data.frame(pairs_data))
}

