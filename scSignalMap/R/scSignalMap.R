##################################################################
##  ______     ______     __  __                                ##
## /\  __ \   /\  ___\   /\ \/\ \                               ##
## \ \  __ \  \ \___  \  \ \ \_\ \                              ##
##  \ \_\ \_\  \/\_____\  \ \_____\                             ##
##   \/_/\/_/   \/_____/   \/_____/                             ##
## @Developed by: Plaisier Lab                                  ##
##   (https://plaisierlab.engineering.asu.edu/)                 ##
##   Arizona State University                                   ##
##   242 ISTB1, 550 E Orange St                                 ##
##   Tempe, AZ  85281                                           ##
## @github: https://github.com/plaisier-lab/scSignalMap         ##
## @Author: Jillian Woodley, Samantha O'Connor, Chris Plaisier  ##
## @License:  GNU GPLv3                                         ##
##                                                              ##
## If this program is used in your analysis please              ##
## mention who built it. Thanks. :-)                            ##
##################################################################


#' Capturing signaling pathways using scRNA-seq data
#'
#' Function to capture signaling pathways using scRNA-seq data.
#'
#'
#' @param seurat_obj: a seurat object must be supplied to classify, no default
#' @param group_by: a meta_data column that will be used to split the cells, no default
#' @param cond_column: meta.data column name for condition information 
#' @param cond_name1: condition column value that compared to cond_name2 for differential expression analysis, must be specified
#' @param avg_log2FC_gte: 
#' @param p_val_adj_lte: 
#' @param species: from which species did the samples originate, either 'human' or 'mouse', defaults to 'human'
#' @param gene_id: what type of gene ID is used, either 'ensembl' or 'symbol', defaults to 'ensembl'
#' @return data.frame with putative interactions
#' @export
map_interactions = function(seurat_obj, group_by, cond_column, cond_name1, avg_log2FC_gte = 0.25, p_val_adj_lte = 0.05, min_pct = 0.1, species='human', gene_id='ensembl') {
    cat('Running scSignalMap:\n')
    

    ## Step 1. Load up the data
    cat('  Loading data...\n')
    # Load MultiNicheNet ligand receptor interactions
    lr_network = read.csv(system.file('extdata', 'lr_network.csv', package='scSignalMap'), header=TRUE)
    lr_pairs = na.omit(lr_network[,c(paste('ligand',species,gene_id,sep='_'),paste('receptor',species,gene_id,sep='_'))])
    #cat(1)
    if(gene_id!='symbol') {
        lig_convert = lr_network[,c(paste('ligand',species,'symbol',sep='_'))]
        names(lig_convert) = lr_network[,c(paste('ligand',species,gene_id,sep='_'))]
        rec_convert = lr_network[,c(paste('receptor',species,'symbol',sep='_'))]
        names(rec_convert) = lr_network[,c(paste('receptor',species,gene_id,sep='_'))]
    }

    # Load up secreted ligands
    secreted = read.csv(system.file('extdata', 'secreted.csv', package='scSignalMap'), header=TRUE)
    colnames(secreted)[1] = 'human_ensembl'
    secreted_ligands = na.omit(secreted[,paste(species,gene_id,sep='_')])


    ## Step 2. Identify marker genes
    cat('  Identifying marker genes...\n')
    
    putative_markers = FindAllMarkers(seurat_obj, group.by=group_by, min.pct=min_pct, verbose=T)
    markers = list()
    for(clust1 in as.character(sort(unique(seurat_obj@meta.data[,group_by])))) {
        markers[[clust1]] = putative_markers %>% dplyr::filter(cluster==clust1 & abs(avg_log2FC)>=avg_log2FC_gte & p_val_adj<=p_val_adj_lte) %>% dplyr::pull(name='gene')
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
    Idents(seurat_obj) = cond_column
    seurat_obj_cond1 = subset(seurat_obj, idents = cond_name1)
    seurat_obj_split = SplitObject(seurat_obj_cond1, split.by=group_by)

    # Build a data.table with statistics
    all_dt = rbindlist(lapply(as.character(sort(unique(seurat_obj_cond1@meta.data[,group_by]))), function(clust1) {
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
    clusts = as.character(sort(unique(seurat_obj_cond1@meta.data[[group_by]])))
    clust_dt = data.table::CJ(Sender = clusts, Receiver = clusts)  # Cartesian product

    # Ligand-Receptor pairs table
    if (gene_id == "symbol") {
      lr_dt = data.table::data.table(
        Ligand   = lr_pairs[, 1],
        Receptor = lr_pairs[, 2]
      )
    } else {
      lr_dt = data.table::data.table(
        Ligand          = lr_pairs[, 1],
        Ligand_Symbol   = lig_convert[lr_pairs[, 1]],
        Receptor        = lr_pairs[, 2],
        Receptor_Symbol = rec_convert[lr_pairs[, 2]]
      )
    }

    # Cartesian product of LR pairs x cluster pairs
    lr_dt[, `:=`(dummy, 1)]
    clust_dt[, `:=`(dummy, 1)]

    pairs_data = lr_dt[clust_dt, on = "dummy", allow.cartesian = TRUE][, `:=`(dummy, NULL)]
    #pairs_data = lr_dt[, cbind(.SD, clust_dt), by = seq_len(nrow(lr_dt))][, `:=`(seq_len, NULL)]

    cat(paste0('    Rows = ',nrow(pairs_data),'\n'))

    cat('  Integrating data...\n')
    steps = 13
    pb = utils::txtProgressBar(min = 0, max = steps, style = 3)
    pairs_data[, `:=`(Ligand_Counts, all_dt[pairs_data, on = .(clust1=Sender, gene=Ligand), counts])]
    utils::setTxtProgressBar(pb, 1)
    pairs_data[, `:=`(Ligand_gte_3, all_dt[pairs_data, on = .(clust1=Sender, gene=Ligand), perc_gte_3])]
    utils::setTxtProgressBar(pb, 2)
    pairs_data[, `:=`(Ligand_gte_10, all_dt[pairs_data, on = .(clust1=Sender, gene=Ligand), perc_gte_10])]
    utils::setTxtProgressBar(pb, 3)
    pairs_data[, `:=`(Ligand_Cells_Exp, all_dt[pairs_data, on = .(clust1=Sender, gene=Ligand), perc_gt_0])]
    utils::setTxtProgressBar(pb, 4)
    pairs_data[, `:=`(Ligand_Avg_Exp, all_dt[pairs_data, on = .(clust1=Sender, gene=Ligand), avg_exp])]
    utils::setTxtProgressBar(pb, 5)
    pairs_data[, `:=`(Ligand_Cluster_Marker, mapply(function(sender, ligand) { ligand %in% markers[[sender]] }, pairs_data$Sender, pairs_data$Ligand))]
    utils::setTxtProgressBar(pb, 6)
    pairs_data[, `:=`(Ligand_secreted, sapply(pairs_data[['Ligand']], function(x) { x %fin% secreted_ligands }))]
    utils::setTxtProgressBar(pb, 7)
    pairs_data[, ':='(Receptor_Counts, all_dt[pairs_data, on = .(clust1=Receiver, gene=Receptor), counts])]
    utils::setTxtProgressBar(pb, 8)
    pairs_data[,`:=`(Receptor_gte_3, all_dt[pairs_data, on = .(clust1=Receiver, gene=Receptor), 'perc_gte_3'])]
    utils::setTxtProgressBar(pb, 9)
    pairs_data[,`:=`(Receptor_gte_10, all_dt[pairs_data, on = .(clust1=Receiver, gene=Receptor), 'perc_gte_10'])]
    utils::setTxtProgressBar(pb, 10)
    pairs_data[,`:=`(Receptor_Cells_Exp, all_dt[pairs_data, on = .(clust1=Receiver, gene=Receptor), perc_gt_0])]
    utils::setTxtProgressBar(pb, 11)
    pairs_data[, `:=`(Receptor_Avg_Exp, all_dt[pairs_data, on = .(clust1=Receiver, gene=Receptor), avg_exp])]
    utils::setTxtProgressBar(pb, 12)
    pairs_data[, `:=`(Receptor_Cluster_Marker, mapply(function(receiver, receptor) { receptor %in% markers[[receiver]] }, pairs_data$Receiver, pairs_data$Receptor))]
    utils::setTxtProgressBar(pb, 13)
    close(pb)
    cat('Done.\n')

    return(pairs_data)
}


##' Perform Differential Expression for a Desired Celltype
#'
#' This function identifies differentially expressed (DE) genes between two conditions 
#' within a specific celltype of interest. It uses Seurat's `FindMarkers` function with 
#' customizable log2FC and adjusted p-value cutoffs. The results are filtered and 
#' annotated with gene symbols from Ensembl IDs. 
#'
#' @param seurat_obj: seurat object containing scRNA data
#' @param prep_SCT: specify whether the PrepSCTFindMarker function needs to be applied, defaults to FALSE
#' @param cond_column: meta.data column name for condition information 
#' @param cond_name1: condition column value that compared to cond_name2 for differential expression analysis, must be specified
#' @param cond_name2: condition column value that compared to cond_name1 for differential expression analysis, must be specified
#' @param celltype_column: meta.data column name for celltype information
#' @param celltype_name: celltype column value of interest, must be specified 
#' @param FC_cutoff: desired cutoff for log2FC values using >= the absolute value, default is 0.3
#' @param adj_p_val_cutoff: desired cutoff for the adjusted p value, default is 0.05
#' @return A data frame of DE genes between conditions for specified cell type, including Ensembl ID, gene symbol, log2FC, and adjusted p-value.
#' @export
find_markers_btwn_cond_for_celltype = function(seurat_obj = NULL, prep_SCT = FALSE, cond_column = NULL, cond_name1 = NULL, cond_name2 = NULL, celltype_column = NULL, celltype_name = NULL, FC_cutoff = 0.3, adj_p_val_cutoff = 0.05, ensdb = 'EnsDb.Hsapiens.v86') {

    message("Preparing to run FindMarkers...")
    if(prep_SCT==TRUE) {
        seurat_obj = PrepSCTFindMarkers(seurat_obj)
    }
    seurat_obj = PrepSCTFindMarkers(seurat_obj)

    message("Subsetting and setting identities...")
    cells.1 = rownames(seurat_obj@meta.data %>% dplyr::filter((!!sym(celltype_column) == celltype_name) & (!!sym(cond_column) == cond_name1)))
    cells.2 = rownames(seurat_obj@meta.data %>% dplyr::filter((!!sym(celltype_column) == celltype_name) & (!!sym(cond_column) == cond_name2)))

    message("Running FindMarkers...")
    de_cells = FindMarkers(seurat_obj, ident.1=cells.1, ident.2=cells.2)

    message("Filtering DE genes by log2FC and adjusted p-value...")
    de_cond_celltype = de_cells %>%
                       dplyr::filter((abs(avg_log2FC) >= FC_cutoff) & (p_val_adj <= adj_p_val_cutoff))

    message("Adding gene symbols...")
    ensembl_ids = rownames(de_cond_celltype)
    
    gene_symbols = AnnotationDbi::select(eval(parse(text=ensdb)), keys = ensembl_ids, keytype = "GENEID", columns = c("SYMBOL"), multiVals='asNA')
    rownames(gene_symbols) = gene_symbols[,'GENEID']

    de_cond_celltype = data.frame(ensembl_id = ensembl_ids,
                                  gene_symbol = gene_symbols[ensembl_ids,'SYMBOL'],
                                  de_cond_celltype,
                                  stringsAsFactors = FALSE)

    de_cond_celltype = de_cond_celltype[, c("ensembl_id", "gene_symbol", setdiff(colnames(de_cond_celltype), c("ensembl_id", "gene_symbol")))]

    return(de_cond_celltype)
}


##' Identify up-regulated receptors
#'
#' This function identifies upregualted receptors from a given DE gene table using a chosen log2FC cutoff
#'
#' @param de_condition_filtered: differentially expressed genes output from find_markers_btwn_cond_for_celltype function
#' @param FC_cutoff: desired cutoff for log2FC values using >= the absolute value, default is 0.3
#' @return A dataframe with identified DE genes and their log2FC from previously chosen condition
#' @export
find_upreg_receptors = function(de_condition_filtered= NULL, FC_cutoff = 0.3, species = 'human') {

    message("Loading ligand-receptor information")
    # Load MultiNicheNet ligand receptor interactions
    lr_network = read.csv(system.file('extdata', 'lr_network.csv', package='scSignalMap'), header=TRUE)
    receptor_ensembl = lr_network[,paste('receptor',species,'ensembl',sep='_')]
    receptor_symbol = lr_network[,paste('receptor',species,'symbol',sep='_')]
    receptor_genes = unique(na.omit(receptor_ensembl))
    ensembl_to_symbol = setNames(receptor_symbol, receptor_ensembl)

    message("Filter for upregulated receptors")
    upreg_receptors = de_condition_filtered %>%
                      dplyr::filter(ensembl_id %in% receptor_genes & avg_log2FC >= FC_cutoff)
    upreg_receptors$gene_symbol = ensembl_to_symbol[upreg_receptors$ensembl_id]

    upreg_receptors = upreg_receptors[, c("gene_symbol", setdiff(names(upreg_receptors), "gene_symbol"))]  

    return(upreg_receptors)
}


##' Identify down-regulated receptors
#'
#' This function identifies down-regualted receptors from a given DE gene table using a chosen log2FC cutoff
#'
#' @param de_condition_filtered: differentially expressed genes output from find_markers_btwn_cond_for_celltype function
#' @param FC_cutoff: desired cutoff for log2FC values using <= the absolute value, default is 0.3
#' @return A dataframe with identified DE genes and their log2FC from previously chosen condition
#' @export
find_downreg_receptors = function(de_condition_filtered= NULL, FC_cutoff = 0.3, species = 'human') {

    message("Loading ligand-receptor information")
    # Load MultiNicheNet ligand receptor interactions
    lr_network = read.csv(system.file('extdata', 'lr_network.csv', package='scSignalMap'), header=TRUE)
    receptor_ensembl = lr_network[,paste('receptor',species,'ensembl',sep='_')]
    receptor_symbol = lr_network[,paste('receptor',species,'symbol',sep='_')]
    receptor_genes = unique(na.omit(receptor_ensembl))
    ensembl_to_symbol = setNames(receptor_symbol, receptor_ensembl)

    message("Filter for downregulated receptors")
    downreg_receptors = de_condition_filtered %>%
                      dplyr::filter(ensembl_id %in% receptor_genes & avg_log2FC <= -FC_cutoff)
    downreg_receptors$gene_symbol = ensembl_to_symbol[downreg_receptors$ensembl_id]

    downreg_receptors = downreg_receptors[, c("gene_symbol", setdiff(names(downreg_receptors), "gene_symbol"))]  

    return(downreg_receptors)
}


##' Filter identified ligand-receptor interactions
#'
#' This function filters interactions from a scSignalMap output based on celltype and ligand secretion
#'
#' @param interactions: ligand-receptor interactions dataframe from MapInteractions function of scSignalMap
#' @param sender_celltypes: celltype of cell(s) that are "senders" in the cell-to-cell communication interaction, ex. "<celltype>" or c("celltype1", "celltype2")
#' @param receiver_celltypes: celltype of cell(s) that are "receivers" in the cell-to-cell communication interaction, ex. "<celltype>" or c("celltype1", "celltype2")  
#' @param secreted_lig: logical statement if the ligand in the interaction is secreted by the sending cell, default is TRUE
#' @return A data frame of filtered interactions
#' @export
filter_lr_interactions = function(interactions = NULL, sender_celltypes = NULL, receiver_celltypes = NULL, secreted_lig = TRUE) {

      message("Filtering scSignalMap interactions")
      
      if (secreted_lig) {
          interactions_filtered = interactions %>%
                                  dplyr::filter(Ligand_secreted == TRUE,
                                  Receiver %in% receiver_celltypes,
                                  Sender %in% sender_celltypes,
                                  Sender != Receiver)
      } else {
          interactions_filtered = interactions %>%
                                  dplyr::filter(Receiver %in% receiver_celltypes,
                                  Sender %in% sender_celltypes,
                                  Sender != Receiver)
      }
      
      return(interactions_filtered)
}


##' Intersect DE receptors with filtered ligand-receptor interactions
#'
#' This function takes previously identified DE receptor genes from chosen condition and identifies overlapping interactions
#'
#' @param de_receptors: name for output file containing DE receptors
#' @param interactions: ligand-receptor interactions dataframe from map_interactions 
#' @return A data frame of idendified DE receptor genes found in previously idendified interactions list
#' @export
intersect_de_receptors_with_lr_interactions = function(de_receptors = NULL, interactions = NULL) {

    message("Intersect DEG receptors with filtered interactions")

    de_receptors_filtered_and_compared = de_receptors %>%
                                            dplyr::filter(gene_symbol %in% interactions$Receptor_Symbol)

return(de_receptors_filtered_and_compared)
}


##' Identify Enriched Pathways from DE Genes Using Enrichr
#'
#' This function performs pathway enrichment analysis on differentially expressed (DE) genes 
#' using the Enrichr. It filters enrichment results by an adjusted p-value cutoff and 
#' saves both individual and combined results for selected databases. The function uses a 
#' user-provided background gene set (from a Seurat object) and compares enrichment results 
#' to input DE genes to highlight relevant pathways.
#'
#' @param seurat_obj: seurat object containing genes for Enrichr to use as background, likely to be seurat object used as input for scSignalMap's MapInteractions function, must provide input
#' @param de_condition_filtered: DE genes output from find_markers_btwn_cond_for_celltype function, used to extract genes for Enrichr use
#' @param enrichr_databases: databases that Enrichr will use for pathway idenification, default includes c("BioCarta_2016", "GO_Biological_Process_2025", "KEGG_2021_Human", "NCI-Nature_2016", "WikiPathways_2024_Human")
#' @param adj_p_val_method: method for p-value adjustments, default is Benjamini-Hochberg (BH) 
#' @param adj_p_val_cutoff: desired cutoff for the adjusted p-value, default is 0.05
#' @return A data frame containing identified pathways, associated statistical values, common genes, etc.
#' @export
find_enriched_pathways = function(seurat_obj = NULL, de_condition_filtered = NULL, enrichr_databases = c("BioCarta_2016", "GO_Biological_Process_2025", "KEGG_2021_Human", "NCI-Nature_2016", "WikiPathways_2024_Human"), adj_p_val_method = "BH", adj_p_val_cutoff = 0.05, ensdb = 'EnsDb.Hsapiens.v86') {

    # Get a list of genes for functional enrichment
    genes = unique(as.character(de_condition_filtered$gene_symbol))

    # Named list to convert gene names
    convert_me = genes
    names(convert_me) = toupper(genes)

    # Generate a background gene list based on all genes in Seurat object
    background_genes = rownames(seurat_obj[["RNA"]])
    background_genes = AnnotationDbi::select(eval(parse(text=ensdb)), keys = background_genes, keytype = "GENEID", columns = c("SYMBOL"))

    # Conduct enrichment analysis using Enrichr
    enrichment_results = enrichR::enrichr(genes, enrichr_databases, background = background_genes[,'SYMBOL'])

    # Add adjusted p-value
    for (db in names(enrichment_results)) {
        data = enrichment_results[[db]]
        data$Adjusted.P.value = p.adjust(data$P.value, method = adj_p_val_method)
        enrichment_results[[db]] = data
    }

    # Make enrichment results directory and save our files for each database tested
    if (!dir.exists("enrichr_results")) {
        dir.create("enrichr_results")
    }
    for (db in names(enrichment_results)) {
        output = paste0("enrichr_results/", db, ".csv")
        write.csv(enrichment_results[[db]], file = output, row.names = FALSE)
    }

    # Combine the enrichr results for all databases
    enrichment_results_combined = bind_rows(
        lapply(names(enrichment_results), function(db) {
            df = enrichment_results[[db]]
            df$database = db
            df })
    )

    # Rename the Genes so they are correct for mouse, and shouldn't hurt human
    enrichr_results = enrichment_results_combined %>% mutate(tmp = paste(convert_me[strsplit(Genes, ";")[[1]]], collapse=';'))

    # Filter down to only the significantly enriched pathways
    enrichr_results = enrichr_results[enrichr_results$Adjusted.P.value < adj_p_val_cutoff, ]

    return(enrichr_results)
}


##' Run Full scSignalMap Pipeline 
#'
#' This function performs a full ligand–receptor signaling analysis  
#' using scSignalMap and Seurat. It performs ligand–receptor interaction mapping, 
#' differential expression analysis, receptor filtering, receptor–ligand 
#' intersection, and pathway enrichment analysis in one workflow.  
#'
#' @param seurat_obj: seurat object containing scRNA-seq data
#' @param prep_SCT: logical; whether to run PrepSCTFindMarkers before differential expression, default = TRUE
#' @param cond_column: meta.data column in Seurat object containing condition labels
#' @param cond_name1: condition column value that compared to cond_name2 for differential expression analysis, must be specified
#' @param cond_name2: condition column value that compared to cond_name1 for differential expression analysis, must be specified
#' @param celltype_column: meta.data column in Seurat object for celltype information
#' @param celltype_name: celltype column value of interest, must be specified
#' @param sender_celltypes: celltype of cell(s) that are "senders" in the cell-to-cell communication interaction, ex. "<celltype>" or c("celltype1", "celltype2")
#' @param receiver_celltypes: celltype of cell(s) that are "receivers" in the cell-to-cell communication interaction, ex. "<celltype>" or c("celltype1", "celltype2")
#' @param secreted_lig: logical; whether to restrict to secreted ligands in filtering, default = TRUE
#' @param FC_cutoff: log2 fold change cutoff for differential expression and receptor filtering, default = 0.3
#' @param adj_p_val_cutoff: adjusted p-value cutoff for differential expression and pathway enrichment, default = 0.05
#' @param enrichr_databases: databases that Enrichr will use for pathway idenification, default includes c("BioCarta_2016", "GO_Biological_Process_2025", "KEGG_2021_Human", "NCI-Nature_2016", "WikiPathways_2024_Human")
#' @param adj_p_val_method: method for p-value adjustments, default = "BH" 
#' @return a list containing ligand–receptor interactions, DE genes, upregulated receptors, filtered interactions, intersected receptors, and pathway enrichment results.
#' @export
run_full_scSignalMap_pipeline = function(seurat_obj = NULL, prep_SCT = TRUE, cond_column = NULL, cond_name1 = NULL, cond_name2 = NULL, celltype_column = NULL, celltype_name = NULL, sender_celltypes = NULL, receiver_celltypes = NULL, secreted_lig = TRUE, FC_cutoff = 0.3, adj_p_val_cutoff = 0.05, enrichr_databases = c("BioCarta_2016", "GO_Biological_Process_2025", "KEGG_2021_Human", "NCI-Nature_2016", "WikiPathways_2024_Human"), adj_p_val_method = "BH", ensdb = 'EnsDb.Hsapiens.v86', species='human') {

  #####################
  ### Run pipeline  ###
  #####################
  message("Running map_interactions...")
  LR_interactions = MapInteractions(seurat_obj, 
                                    group_by = celltype_column,
                                    cond_column = cond_column,
                                    cond_name1= cond_name1,
                                    species=species)

  message("Finding DE genes...")
  de_cond_celltype = find_markers_btwn_cond_for_celltype(
      seurat_obj = seurat_obj,
      prep_SCT = prep_SCT,
      cond_column = cond_column,
      cond_name1 = cond_name1,
      cond_name2 = cond_name2,
      celltype_column = celltype_column,
      celltype_name = celltype_name,
      FC_cutoff = FC_cutoff,
      adj_p_val_cutoff = adj_p_val_cutoff,
      ensdb = ensdb)

  message("Finding upregulated receptors...")
  upreg_receptors = find_upreg_receptors(
      de_condition_filtered = de_cond_celltype,
      FC_cutoff = FC_cutoff, species=species)

  message("Finding upregulated receptors...")
  downreg_receptors = find_downreg_receptors(
      de_condition_filtered = de_cond_celltype,
      FC_cutoff = FC_cutoff, species=species)

  message("Filtering LR interactions...")
  interactions_filtered = filter_lr_interactions(
      interactions = LR_interactions,
      sender_celltypes = sender_celltypes,
      receiver_celltypes = receiver_celltypes,
      secreted_lig = secreted_lig)

  message("Intersecting receptors with interactions (up)...")
  upreg_receptors_filtered_and_compared = intersect_de_receptors_with_lr_interactions(
      upreg_receptors = upreg_receptors,
      interactions = interactions_filtered)

  message("Intersecting receptors with interactions (down)...")
  downreg_receptors_filtered_and_compared = intersect_de_receptors_with_lr_interactions(
      downreg_receptors = downreg_receptors,
      interactions = interactions_filtered)

  message("Integrate up and down regualted receptors...")
  combined_receptors_filtered_and_compared = rbind(upreg_receptors_filtered_and_compared, downreg_receptors_filtered_and_compared)
  combined_receptors_filtered_and_compared[,'Feeback'] = ifelse(combined_receptors_filtered_and_compared[,'avg_log2FC']>0, 'Amplification','Adaptation')

  message("Running pathway enrichment...")
  enrichr_results = find_enriched_pathways(
      seurat_obj = seurat_obj,
      de_condition_filtered = de_cond_celltype,
      enrichr_databases = enrichr_databases,
      adj_p_val_method = adj_p_val_method,
      adj_p_val_cutoff = adj_p_val_cutoff,
      ensdb = ensdb)

  ###################################
  ### Return all results together ###
  ###################################
  return(list(
      LR_interactions = LR_interactions,
      de_cond_celltype = de_cond_celltype,
      upreg_receptors = upreg_receptors,
      downreg_receptors = downreg_receptors,
      interactions_filtered = interactions_filtered,
      upreg_receptors_filtered_and_compared = upreg_receptors_filtered_and_compared,
      downreg_receptors_filtered_and_compared = downreg_receptors_filtered_and_compared,
      combined_receptors_filtered_and_compared = combined_receptors_filtered_and_compared,
      enrichr_results = enrichr_results))
}


##' Create Master Interaction List
#'
#' This function creates master interaction list by combining DE ligands/receptors, Enrichr results, and scSignalMap interactions
#'
#' @param enrichr_results Data frame of Enrichr pathway enrichment results, default is results$enrichr_results 
#' @param de_receptors Data frame of up- and down-regulated receptors, default is results$combined_receptors_filtered_and_compared
#' @param scSignalMap_data_filtered Data frame of filtered interactions, default is results$interactions_filtered
#' @return A data frame containing merged ligand/receptor info, enrichment, and interaction data
#' @export
create_master_interaction_list = function(
  enrichr_results = results$enrichr_results,
  de_receptors = results$combined_receptors_filtered_and_compared,
  scSignalMap_data_filtered = results$interactions_filtered) {

    ## Step 1: Clean Enrichr results
    enrichr_results = enrichr_results[, !(names(enrichr_results) %in% c("Old.P.value", "Old.Adjusted.P.value"))]

    # Break appart the genes
    genes = sapply(enrichr_results$Genes, function(x) { strsplit(x,';')[[1]] })
    names(genes) = enrichr_results$Term

    # Filtering down to terms
    matched = list()
    for(term1 in names(genes)) {
        intersect1 = intersect(genes[[term1]],unique(de_receptors$gene_symbol))
        if(length(intersect1)>0) {
            matched[[term1]] = intersect1
        }
    }

    # Make the master list
    master_list = vector("list", length = 0)
    for (term1 in names(matched)) {
        cur_term_df = enrichr_results %>%
          dplyr::filter(Term == term1)
        for (rec1 in matched[[term1]]) {
            rec1_df = scSignalMap_data_filtered %>%
              dplyr::filter(Receptor_Symbol == rec1)
            deg_res = de_receptors %>%
              dplyr::filter(gene_symbol == rec1)
            if (nrow(rec1_df) > 0) {
                combined_df = bind_cols(cur_term_df[rep(1, nrow(rec1_df)), ], rec1_df, deg_res[,c('p_val','p_val_adj','avg_log2FC','Feedback','pct.1','pct.2')])
                master_list[[length(master_list) + 1]] = combined_df
            }
        }
    }

    master_list = bind_rows(master_list)

    return(master_list)
}
