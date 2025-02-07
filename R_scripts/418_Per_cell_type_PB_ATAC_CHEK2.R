.libPaths()
.libPaths(new = c("/home/manuel.tardaguila/conda_envs/multiome_NEW_downstream_analysis/lib/R/library"))
.libPaths()
# sessionInfo()
Sys.setenv(RETICULATE_PYTHON="/home/manuel.tardaguila/conda_envs/multiome_NEW_downstream_analysis/bin/python")

suppressMessages(library("optparse"))
suppressMessages(library(Seurat))
suppressMessages(library(Signac))
suppressMessages(library(dplyr)) 
suppressMessages(library(ggplot2)) 
suppressMessages(library(Matrix)) 
suppressMessages(library(data.table)) 
suppressMessages(library(ggpubr)) 
suppressMessages(library(ggplot2))
suppressMessages(library(pheatmap))
suppressMessages(library("qlcMatrix"))
suppressMessages(library("cowplot"))
suppressMessages(library("RColorBrewer"))
suppressMessages(library("plyr"))
suppressMessages(library("forcats"))
suppressMessages(library('ggeasy'))
suppressMessages(library('dplyr'))
suppressMessages(library("svglite"))
suppressMessages(library("ape"))
suppressMessages(library("ggforce"))
suppressMessages(library("tidyr"))
suppressMessages(library("edgeR"))
suppressMessages(library("apeglm"))
suppressMessages(library("DESeq2"))
suppressMessages(library("tibble")) 


opt = NULL

options(warn = 0)


data_wrangling = function(option_list)
{
  
  opt_in = option_list
  opt <<- option_list
  
  cat("All options:\n")
  printList(opt)
  
  
  #### READ and transform type ----
  
  type = opt$type
  
  cat("TYPE_\n")
  cat(sprintf(as.character(type)))
  cat("\n")
  
  #### READ and transform out ----
  
  out = opt$out
  
  cat("out_\n")
  cat(sprintf(as.character(out)))
  cat("\n")
  
  #### READ and transform out ----
  
  cell_type_sel = opt$cell_type_sel
  
  cat("cell_type_sel_\n")
  cat(sprintf(as.character(cell_type_sel)))
  cat("\n")


  
 
  #### READ adata ----
  
  adata<-readRDS(file=opt$adata)
  

  # cat("Seurat_object_0\n")
  # cat(str(adata))
  # cat("\n")
  
  # Create a genotype factor ----------------------
  
  adata@meta.data$Assigned_GFPgenotype_integral<-NA
  
  adata@meta.data$Assigned_GFPgenotype_integral[which(adata@meta.data$Assigned_GFPbc_integral%in%c('chrGFP_WTA','chrGFP_WTB','chrGFP_WTC'))]<-'wt'
  adata@meta.data$Assigned_GFPgenotype_integral[which(adata@meta.data$Assigned_GFPbc_integral%in%c('chrGFP_rs1','chrGFP_rs2','chrGFP_rs3'))]<-'CHEK2 T/T'
  adata@meta.data$Assigned_GFPgenotype_integral[which(adata@meta.data$Assigned_GFPbc_integral%in%c('chrGFP_R882H1','chrGFP_R882H2','chrGFP_R882H3'))]<-'DNMT3A R882H'
  adata@meta.data$Assigned_GFPgenotype_integral[which(adata@meta.data$Assigned_GFPbc_integral%in%c('chrGFP_rs_R882H1','chrGFP_rs_R882H2','chrGFP_rs_R882H3'))]<-'Double mutants'
  
  

  
  adata@meta.data$Assigned_GFPgenotype_integral<-factor(adata@meta.data$Assigned_GFPgenotype_integral,
                                                        levels=c('wt','CHEK2 T/T','DNMT3A R882H','Double mutants'),
                                                        ordered=T)
  
  # Subset my object for cells with genotype assigned ---------------------
  
  
  adata_geno<-subset(adata, Assigned_GFPgenotype_integral != "NA")

  # Subset for cell type ---------------------------------

  ct<-as.character(unique(adata_geno@meta.data$current_anot[which(as.numeric(adata_geno@meta.data$current_anot) == cell_type_sel)]))

  adata_geno_ct_sel<-subset(adata_geno, current_anot == ct)

  
  ## Create a variable to correct for clone line -----------------------------------
  
  adata_geno_ct_sel$sample_id<-adata_geno_ct_sel$Assigned_GFPbc_integral
  
  cat(sprintf(as.character(names(summary(adata_geno_ct_sel$sample_id)))))
  cat("\n")
  cat(sprintf(as.character(summary(adata_geno_ct_sel$sample_id))))
  cat("\n")
  
  ## Extract ATAC counts not normalized -----------------------
  
  
  matrix_ATAC<-GetAssayData(object = adata_geno_ct_sel, assay = "ATAC", layer = "counts")
  
  cat("matrix_ATAC\n")
  cat(str(matrix_ATAC))
  cat("\n")
  
  ## Extract metadata -------------------------------------------
  
  metadata<-droplevels(adata_geno_ct_sel[[]])
  
  cat(sprintf(as.character(names(summary(metadata$time_point)))))
  cat("\n")
  cat(sprintf(as.character(summary(metadata$time_point))))
  cat("\n")
  
  
  cat(sprintf(as.character(names(summary(metadata$current_anot)))))
  cat("\n")
  cat(sprintf(as.character(summary(metadata$current_anot))))
  cat("\n")
  
  
  cat(sprintf(as.character(names(summary(as.factor(metadata$Assigned_GFPbc_integral))))))
  cat("\n")
  cat(sprintf(as.character(summary(as.factor(metadata$Assigned_GFPbc_integral)))))
  cat("\n")
  
  cat(sprintf(as.character(names(summary(as.factor(metadata$sample_id))))))
  cat("\n")
  cat(sprintf(as.character(summary(as.factor(metadata$sample_id)))))
  cat("\n")
  
  
  ## Create a new Seurat object with the ATAC and the metadata -------------------
  
  ATAC_object <- CreateSeuratObject(counts = matrix_ATAC, assay = "RNA",
                                    meta.data=metadata)
  
  ## Aggregate by sample_id (clone line) and current_anot (cell type) ----------------------------------
  
  
  cell_type_names <- levels(metadata[,which(colnames(metadata) == 'current_anot')])
  
  cat("cell_type_names_0\n")
  cat(str(cell_type_names))
  cat("\n")
  
  sample_names <- levels(metadata[,which(colnames(metadata) == 'sample_id')])
  
  cat("sample_names_0\n")
  cat(str(sample_names))
  cat("\n")
  
  groups <- metadata[,c(which(colnames(metadata) == 'sample_id'),which(colnames(metadata) == 'current_anot'))]
  
  cat("groups_0\n")
  cat(str(groups))
  cat("\n")
  
  aggr_counts <- Seurat2PB(ATAC_object, sample="sample_id", cluster="current_anot")
  
  cat("aggr_counts_0\n")
  cat(str(aggr_counts))
  cat("\n")
  
  
  ## Go to a list format for counts --------------------------------------------------------
  
  
  ## Initiate empty list
  counts_ls <- list()
  
  DEBUG<-0
  
  
  for (i in 1:length(cell_type_names)) {
    
    cell_type_names[i]
    
    ## Extract indexes of columns in the global matrix that match a given cluster
    column_idx <- which(tstrsplit(colnames(aggr_counts), "_cluster")[[2]] == cell_type_names[i])
    
    sub_aggr<- aggr_counts[, column_idx]
    
    if(DEBUG == 1)
    {
      cat("sub_aggr_0\n")
      cat(str(sub_aggr))
      cat("\n")
      
    }
    
    ## Store corresponding sub-matrix as one element of a list
    counts_ls[[i]] <-sub_aggr
    names(counts_ls)[i] <- cell_type_names[i]
    
    #break
    
  }
  
  
  ## Create a metadata list -------------------------------------------------------
  
  #### Create group level variables -----
  
  # Extract sample-level variables
  metadata_NEW <- metadata %>% 
    as.data.frame() %>% 
    dplyr::select(Assigned_GFPgenotype_integral, Assigned_GFPbc_integral, sample_id)
  
  cat("metadata_NEW_0\n")
  cat(str(metadata_NEW))
  cat("\n")
  
  # Exclude duplicated rows
  metadata_NEW <- metadata_NEW[!duplicated(metadata_NEW), ]
  
  
  cat("metadata_NEW_0.5\n")
  cat(str(metadata_NEW))
  cat("\n")
  
  
  # Rename rows
  rownames(metadata_NEW) <- metadata_NEW$sample_id
  
  cat("metadata_NEW_1\n")
  cat(str(metadata_NEW))
  cat("\n")
  
  t <- table(metadata$sample_id,
             metadata$current_anot)
  
  cat("t_0\n")
  cat(str(t))
  cat("\n")
  
  
  ##### Creating metadata list ----------------
  
  ## Initiate empty list
  metadata_ls <- list()
  
  DEBUG<-0
  
  for (i in 1:length(counts_ls)) {
    
    ## Initiate a data frame for cluster i with one row per sample (matching column names in the counts matrix)
    df <- data.frame(current_anot_sample_id = colnames(counts_ls[[i]]))
    
    if(DEBUG == 1){
      
      cat("df_0\n")
      cat(str(df))
      cat("\n")
    }
    
    ## Use tstrsplit() to separate cluster (cell type) and sample IDs
    df$current_anot_id <- tstrsplit(df$current_anot_sample_id, "_cluster")[[2]]
    df$sample_id  <- tstrsplit(df$current_anot_sample_id, "_cluster")[[1]]
    
    
    if(DEBUG == 1){
      
      cat("df_1\n")
      cat(str(df))
      cat("\n")
    }
    
    
    ## Retrieve cell count information for this cluster from global cell count table
    idx <- which(colnames(t) == unique(df$current_anot_id))
    
    if(DEBUG == 1){
      
      cat("idx_0\n")
      cat(str(idx))
      cat("\n")          
    }
    
    cell_counts <- t[, idx]
    
    if(DEBUG == 1){
      
      cat("cell_counts_0\n")
      cat(str(cell_counts))
      cat("\n")          
    }
    
    
    
    
    
    
    ## Remove samples with zero cell contributing to the cluster
    cell_counts <- cell_counts[cell_counts > 0]
    
    if(DEBUG == 1){
      
      cat("cell_counts_1\n")
      cat(str(cell_counts))
      cat("\n")          
    }
    
    
    ## Match order of cell_counts and sample_ids
    sample_order <- match(df$sample_id, names(cell_counts))
    cell_counts <- cell_counts[sample_order]
    
    if(DEBUG == 1){
      
      cat("cell_counts_2\n")
      cat(str(cell_counts))
      cat("\n")          
    }
    
    
    
    ## Append cell_counts to data frame
    df$cell_count <- cell_counts
    
    if(DEBUG == 1){
      
      cat("df_2\n")
      cat(str(df))
      cat("\n")
    }
    
    
    ## Join data frame (capturing metadata specific to cluster) to generic metadata
    df <- plyr::join(df, metadata_NEW, 
                     by = intersect(names(df), names(metadata_NEW)))
    
    if(DEBUG == 1){
      
      cat("df_3\n")
      cat(str(df))
      cat("\n")
    }
    
    ## Update rownames of metadata to match colnames of count matrix, as needed later for DE
    rownames(df) <- df$current_anot_sample_id
    
    if(DEBUG == 1){
      
      cat("df_4\n")
      cat(str(df))
      cat("\n")
    }
    
    ## Store complete metadata for cluster i in list
    metadata_ls[[i]] <- df
    names(metadata_ls)[i] <- unique(df$cluster_id)
    
    # break
    
  }
  
  
  # DESeq2 for DA analysis ----------------------------------------------------------------------------------------------------

  check_1<-all(names(counts_ls) == names(metadata_ls))
  
  cat("check_1\n")
  cat(sprintf(as.character(check_1)))
  cat("\n")

  array_cell_annotations<-names(counts_ls)

DEBUG<-0

Global_results<-data.frame()

Global_normalized<-data.frame()

for(i in 1:length(array_cell_annotations)){

    Results_per_comparison_within_cell_type<-data.frame()


    cell_type_DEF<-array_cell_annotations[i]

    cat("----------------------------------->\t")
    cat(sprintf(as.character(cell_type_DEF)))
    cat("\n")

    idx.sel<-which(names(counts_ls) == cell_type_DEF)
      cluster_counts <- counts_ls[[idx.sel]]

    if(DEBUG == 1){
          cat("cluster_counts_0\n")
          cat(str(cluster_counts))
          cat("\n")
    }

     cell_type_metadata <- metadata_ls[[idx.sel]]

      if(DEBUG == 1){
          cat("cell_type_metadata_0\n")
          cat(str(cell_type_metadata))
          cat("\n")
    }

    
      cell_type_metadata$Genotype<-as.character(cell_type_metadata$Assigned_GFPgenotype_integral)
      cell_type_metadata$Genotype<-factor(cell_type_metadata$Genotype)
      cell_type_metadata$Genotype<-relevel(cell_type_metadata$Genotype, ref='wt') ### Nothing works with ordered factors


       if(DEBUG == 1){
          cat("cell_type_metadata_1\n")
          cat(str(cell_type_metadata))
          cat("\n")
          cat(str(unique(cell_type_metadata$current_anot_id)))
          cat("\n")
    }

    check_2<-all(colnames(cluster_counts) == rownames(cell_type_metadata))

    if(DEBUG == 1){
      cat("check_2\n")
      cat(sprintf(as.character(check_2)))
      cat("\n")
        cat(sprintf(as.character(rownames(cell_type_metadata))))
      cat("\n")
    }


    
    ### DESeq data set with all the columns the contrast will be every column against the reference wt level set before

     dds_NEW <- DESeqDataSetFromMatrix(cluster_counts, 
                                  colData = cell_type_metadata, 
                                  design = ~ Genotype)

    ### Normalize ----------------------------------

    dds_NEW<-estimateSizeFactors(dds_NEW)

    #### Run DESeq LRT comparing to the reduced model

    dds_NEW_lrt <- DESeq(dds_NEW, test = "LRT", reduced = ~ 1)

   
     ###### Obtain all the possible contrasts --------------------------

    possible_contrasts<-colnames(dds_NEW_lrt@modelMatrix)[-1] # -1 because 1 is the Intercept term

    for(iteration_contrasts in 1:length(possible_contrasts)){

        contrast_sel<-possible_contrasts[iteration_contrasts]

        cat("------------------------------------->\t")
        cat(sprintf(as.character(contrast_sel)))
        cat("\n")
        

        ############ obtain the result for this contrast with its specific pvalue coming from Wald test see ?results

        tmp_results<-results(dds_NEW_lrt, test= "Wald", name=contrast_sel)

        #### expand LogFC #########################################
    
        tmp_results <- lfcShrink(dds_NEW_lrt, 
                                     coef = contrast_sel,
                                     res=tmp_results,
                                     type = "apeglm")

          if(DEBUG == 1){
          cat("tmp_results_0\n")
          cat(str(tmp_results))
          cat("\n")
            }


          #### obtain data frame #########################################

          tmp_tb <- as.data.frame(tmp_results %>%
                      data.frame() %>%
                      rownames_to_column(var = "Peak_ID") %>%
                      as_tibble() %>%
                      arrange(padj), stringsAsFactors=F)
                    
                    

        if(DEBUG == 1){
          cat("tmp_tb_0\n")
          cat(str(tmp_tb))
          cat("\n")
            }

        tmp_tb$contrast<-contrast_sel


        if(DEBUG == 1){
          cat("tmp_tb_1\n")
          cat(str(tmp_tb))
          cat("\n")
            }


        Results_per_comparison_within_cell_type<-rbind(tmp_tb,Results_per_comparison_within_cell_type)
    

        

        }#iteration_contrasts in 1:length(possible_contrasts)


    Results_per_comparison_within_cell_type$cell_type<-cell_type_DEF

    Global_results<-rbind(Results_per_comparison_within_cell_type,Global_results)

    #### Extract normalized counts to plot heatmaps later ----------------------------------------

    nor_counts<-as.data.frame(counts(dds_NEW, normalized=TRUE))

    nor_counts$Peak_ID<-row.names(nor_counts)
    row.names(nor_counts)<-NULL

    nor_counts.m<-melt(nor_counts, id.vars='Peak_ID', variable.name="sample_id_current_annot_id", value.name="count")
    nor_counts.m$GFPbc_integral<-gsub("_cluster.+$","",nor_counts.m$sample_id_current_annot_id)
    nor_counts.m$current_anot<-gsub("^.+_cluster","",nor_counts.m$sample_id_current_annot_id)    
    nor_counts.m<-nor_counts.m[,-which(colnames(nor_counts.m) == 'sample_id_current_annot_id')]

    if(DEBUG == 1){
       cat("nor_counts.m_0\n")
         str(nor_counts.m)
         cat("\n")
        }
    

    Global_normalized<-rbind(nor_counts.m, Global_normalized)

    setwd(out)


    write.table(Global_results, file=paste("DA_results_without_time_as_a_covariate_NEW_METHOD_only_linked_peaks_to_DE_genes_",cell_type_DEF,".tsv", sep=''), sep="\t", quote=F, row.names=F)

    write.table(Global_normalized, file=paste("normalized_Peak_counts_NEW_METHOD_only_linked_peaks_to_DE_genes",cell_type_DEF,".tsv",sep=''), sep="\t", quote=F, row.names=F)

    
    
}# i in 1:length(array_cell_annotations
  

}


printList = function(l, prefix = "    ") {
  list.df = data.frame(val_name = names(l), value = as.character(l))
  list_strs = apply(list.df, MARGIN = 1, FUN = function(x) { paste(x, collapse = " = ")})
  cat(paste(paste(paste0(prefix, list_strs), collapse = "\n"), "\n"))
}



#### main script ----

main = function() {
  cmd_line = commandArgs()
  cat("Command line:\n")
  cat(paste(gsub("--file=", "", cmd_line[4], fixed=T),
            paste(cmd_line[6:length(cmd_line)], collapse = " "),
            "\n\n"))
  option_list <- list(
    make_option(c("--adata"), type="character", default=NULL, 
                metavar="type", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--cell_type_sel"), type="numeric", default=NULL, 
                metavar="type", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--type"), type="character", default=NULL, 
                metavar="type", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--out"), type="character", default=NULL, 
                metavar="type", 
                help="Path to tab-separated input file listing regions to analyze. Required.")
  )
  parser = OptionParser(usage = "140__Rscript_v106.R
                        --subset type
                        --TranscriptEXP FILE.txt
                        --cadd FILE.txt
                        --ncboost FILE.txt
                        --type type
                        --out filename",
                        option_list = option_list)
  opt <<- parse_args(parser)
  
  data_wrangling(opt)
 

}


###########################################################################

system.time( main() )
  