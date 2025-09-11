#----docker 
#docker run -it -v "/home/jswoodl2/MN_data:/files" --rm cplaisier/ccafv2_s5_extra

#Install
remotes::install_github('plaisier-lab/scSignalMap/scSignalMap')
install.packages("enrichR")
install.packages("fastmatch")
library(scSignalMap)
library(enrichR)
library(Seurat)
library(dplyr)
library(tidyr)
library(stringr)
library(org.Hs.eg.db)
library(AnnotationDbi)
library(fastmatch)

#_______________________________________________________________


#################################################
### Perform DE for Desired Celltype & Cutoffs ###
#################################################

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
find_markers_btwn_cond_for_celltype = function(
    seurat_obj = NULL,
    prep_SCT = FALSE,
    cond_column = NULL,
    cond_name1 = NULL,
    cond_name2 = NULL,
    celltype_column = NULL,
    celltype_name = NULL,
    FC_cutoff = 0.3,
    adj_p_val_cutoff = 0.05) {

    # Subset cells and set identities for cond2 and cond1
    message("Subsetting and setting identities...")
    cells_to_keep = rownames(seurat_obj@meta.data)[seurat_obj@meta.data[,celltype_column] == celltype_name]
    subset_cells = subset(seurat_obj, cells = cells_to_keep)
    Idents(subset_cells) = subset_cells@meta.data[,cond_column]

    # FindMarkers
    message("Preparing for and running FindMarkers...")
    if(prep_SCT==TRUE) {
        subset_cells = PrepSCTFindMarkers(subset_cells)
    }
    de_cells = FindMarkers(subset_cells, ident.1 = cond_name1, ident.2 = cond_name2)

    # Filter based on log2FC and p-value
    message("Filtering DE genes by log2FC and adjusted p-value...")
    de_cond_celltype = de_cells %>%
    filter((avg_log2FC >= FC_cutoff | avg_log2FC <= -FC_cutoff) & p_val_adj <= adj_p_val_cutoff)

    # Extract Ensembl IDs from rownames and map to gene symbols
    message("Adding gene symbols...")
    ensembl_ids = rownames(de_cond_celltype)
    gene_symbols = mapIds(org.Hs.eg.db,
                          keys = ensembl_ids,
                          column = "SYMBOL",
                          keytype = "ENSEMBL",
                          multiVals = "first")

    # Clean up columns
    de_cond_celltype = data.frame(ensembl_id = ensembl_ids,
                                  gene_symbol = gene_symbols,
                                  de_cond_celltype,
                                  stringsAsFactors = FALSE)
    
    # Reorder to ensembl_id and gene_symbol are first and second columns
    de_cond_celltype = de_cond_celltype[, c("ensembl_id", "gene_symbol", setdiff(colnames(de_cond_celltype), c("ensembl_id", "gene_symbol")))]
    
    # Return results
    return(de_cond_celltype)
}


##################################
### Find Upregualted Receptors ###
##################################

##' Identify upregulated receptors
#'
#' This function identifies upregualted receptors from a given DE gene table using a chosen log2FC cutoff
#'
#' @param de_condition_filtered: differentially expressed genes output from find_markers_btwn_cond_for_celltype function
#' @param FC_cutoff: desired cutoff for log2FC values using >= the absolute value, default is 0.3
#' @return A dataframe with identified DE genes and their log2FC from previously chosen condition
#' @export
find_upreg_receptors = function(
    de_condition_filtered= NULL,
    FC_cutoff = 0.3) {
    
    # Load ligand-receptor list from scSignalMap
    message("Loading ligand-receptor information")
    lr_list = read.csv(system.file('extdata', 'lr_network.csv', package='scSignalMap'))
    receptor_genes = unique(na.omit(lr_list$receptor_human_ensembl))
    ensembl_to_symbol = setNames(lr_list$receptor_human_symbol, lr_list$receptor_human_ensembl)

    # Filter for just upregulated receptors, add gene symbols, organize columns
    message("Filter for upregulated receptors")
    upreg_receptors = de_condition_filtered %>%
                          filter(ensembl_id %in% receptor_genes & avg_log2FC >= FC_cutoff)
    upreg_receptors$gene_symbol = ensembl_to_symbol[upreg_receptors$ensembl_id]

    # Organize the columns to place the gene_symbol first
    upreg_receptors = upreg_receptors[, c("gene_symbol", setdiff(names(upreg_receptors), "gene_symbol"))]  

    # Return results
    return(upreg_receptors)
}


#########################################
## Filter Ligand-Receptor Interactions ##
#########################################


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
filter_lr_interactions = function(
    interactions = NULL,
    sender_celltypes = NULL,
    receiver_celltypes = NULL,
    secreted_lig = TRUE) {
  
    message("Filtering scSignalMap interactions")
    
    if (secreted_lig) {
        interactions_filtered = interactions %>%
            filter(Ligand_secreted == TRUE,
                   Receiver %in% receiver_celltypes,
                   Sender %in% sender_celltypes,
                   Sender != Receiver)
    } else {
        interactions_filtered = interactions %>%
            filter(Receiver %in% receiver_celltypes,
                   Sender %in% sender_celltypes,
                   Sender != Receiver)
    }
    
    return(interactions_filtered)
}



#######################################################################
## Intersect Upregulated Receptors with Ligand-Receptor Interactions ##
#######################################################################

##' Intersect upregulated receptors with filtered ligand-receptor interactions
#'
#' This function takes previously identified DE receptor genes from chosen condition and identifies overlapping interactions
#'
#' @param upreg_receptors: name for output file containing DE receptors
#' @param interactions: ligand-receptor interactions dataframe from map_interactions 
#' @return A data frame of idendified DE receptor genes found in previously idendified interactions list
#' @export
intersect_upreg_receptors_with_lr_interactions = function(
    upreg_receptors = NULL,
    interactions = NULL) {

    message("Intersect upregulated receptors with filtered interactions")
    
    upreg_receptors_filtered_and_compared = upreg_receptors %>%
        filter(gene_symbol %in% interactions$Receptor_Symbol)

    return(upreg_receptors_filtered_and_compared)
}


#_______________________________________________________________


##########################################################
### Pathways Enriched with DE Genes in Receiving Cells ###
##########################################################

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
find_enriched_pathways = function(
    seurat_obj = NULL,
    de_condition_filtered = NULL,
    enrichr_databases = c("BioCarta_2016", "GO_Biological_Process_2025", "KEGG_2021_Human", "NCI-Nature_2016", "WikiPathways_2024_Human"),
    adj_p_val_method = "BH",
    adj_p_val_cutoff = 0.05) {

    # Extract gene symbols to input in Enrichr
    genes = unique(as.character(de_condition_filtered$gene_symbol))

    # Load in background genes
    background_genes = rownames(seurat_obj[["RNA"]])
    background_genes = mapIds(org.Hs.eg.db, 
                              keys = background_genes, 
                              column = "SYMBOL", 
                              keytype = "ENSEMBL", 
                              multiVals = "first")

    # Run Enrichr for both up and down regulated genes:
    enrichment_results = enrichr(genes, enrichr_databases, background = background_genes)

    for (db in names(enrichment_results)) {
        data = enrichment_results[[db]]

        # Add P.value with adjusted p-values 
        data$Adjusted.P.value = p.adjust(data$P.value, method = adj_p_val_method)

        enrichment_results[[db]] = data
    }

    # Save each db's results as its own .csv
    for (db in names(enrichment_results)) {
        output = paste0("enrichr_results/", db, ".csv")
        write.csv(enrichment_results[[db]], file = output, row.names = FALSE)
    }

    # Combine results:
    enrichment_results_combined = bind_rows(
        lapply(names(enrichment_results), function(db) {
            df = enrichment_results[[db]]
            df$database = db
            df })
    )

    ## Final Comparison:

    # Load in and extract de genes to compare enrichr results to:
    de_genes = unique(de_condition_filtered$gene_symbol)

    # Keep enrichr results that contain at least one de gene in the "Genes" column
    enrichr_results = filter(enrichment_results_combined, grepl(paste(de_genes, collapse="|"), Genes))

    # Keep rows with adj.p < 0.05
    enrichr_results = enrichr_results[enrichr_results$Adjusted.P.value < adj_p_val_cutoff, ]

    # Return results
    return(enrichr_results)
}


###################################################################
############# ---------- Test Out Analysis ---------- #############
###################################################################

# Set working directory
setwd('/files')


# Load seurat object
seurat_obj = readRDS('seurat_objects/MN_big_int.rds')

# If cell types are wanted instead of numbered clusters use:
# and change group_by in call below from 'seurat_clusters' to 'celltype'
seurat_obj$celltype = Idents(seurat_obj)


## Run function and save as .csv (add gene_id='symbol' if needed)
LR_interactions = MapInteractions(seurat_obj,
                                  'celltype')

write.csv(LR_interactions, "MN_big_int_scSignalMap_interactions.csv")


## Find DE genes between MN2 and MN1 for GB3_Tumor cells
de_cond_celltype = find_markers_btwn_cond_for_celltype(
    seurat_obj = seurat_obj,
    prep_SCT = TRUE,
    cond_column = "orig.ident",
    cond_name1 = "MN2",
    cond_name2 = "MN1",
    celltype_column = "celltype",
    celltype_name = "GB3_Tumor",
    FC_cutoff = 0.3,
    adj_p_val_cutoff = 0.05)

write.csv(de_cond_celltype, "de_btwn_MN2_vs_MN1_GB3_Tumor.csv", row.names = FALSE)


## Find upregulated receptors
upreg_receptors = find_upreg_receptors(
    de_condition_filtered = de_cond_celltype,
    FC_cutoff = 0.3)

write.csv(upreg_receptors, "upregulated_receptors_btwn_MN2_vs_MN1_GB3_Tumor.csv", row.names = FALSE)

## Filter LR Interactions
interactions_filtered = filter_lr_interactions(
    interactions = LR_interactions,
    sender_celltypes = c("Astrocyte", "HUVEC"),
    receiver_celltypes = "GB3_Tumor",
    secreted_lig = TRUE)


## Intersect upregulated receptors with ligand-receptor pairs
upreg_receptors_filtered_and_compared = intersect_upreg_receptors_with_lr_interactions(
    upreg_receptors = upreg_receptors,
    interactions = interactions_filtered
)
write.csv(upreg_receptors_filtered_and_compared, "upreg_receptors_filtered_and_compared.csv", row.names = FALSE)

## Find pathways that are enriched with DE genes that include the upregulated receptor
enrichr_results = find_enriched_pathways(
    seurat_obj = seurat_obj,
    de_condition_filtered = de_cond_celltype,
    enrichr_databases = c("BioCarta_2016", "GO_Biological_Process_2025", "KEGG_2021_Human", "NCI-Nature_2016", "WikiPathways_2024_Human"),
    adj_p_val_method = "BH",
    adj_p_val_cutoff = 0.05)

write.csv(enrichr_results, "enrichr_results.csv", row.names = FALSE)


