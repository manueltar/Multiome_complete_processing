
.libPaths()
.libPaths(new = c("/home/manuel.tardaguila/conda_envs/multiome_QC_DEF/lib/R/library"))
.libPaths()
# sessionInfo()

Sys.setenv(RETICULATE_PYTHON="/home/manuel.tardaguila/conda_envs/multiome_QC_DEF/bin/python")
library(reticulate)
reticulate::use_python("/home/manuel.tardaguila/conda_envs/multiome_QC_DEF/bin/python")
reticulate::use_condaenv("/home/manuel.tardaguila/conda_envs/multiome_QC_DEF")
reticulate::py_module_available(module='leidenalg')
reticulate::import('leidenalg')
suppressMessages(library("optparse"))
suppressMessages(library(hdf5r))
suppressMessages(library(Seurat))
suppressMessages(library(Signac))
suppressMessages(library(dplyr))
suppressMessages(library(ggplot2))
suppressMessages(library(Matrix))
suppressMessages(library(data.table))
suppressMessages(library(ggpubr))
suppressMessages(library(ggplot2))
suppressMessages(library(scDblFinder))
suppressMessages(library("tidyr"))
suppressMessages(library("tibble"))
# suppressMessages(library("biovizBase"))
# suppressMessages(library("patchwork"))
suppressMessages(library(glmGamPoi))


refine_metadata_levels <- function(seurat_data){
  for (i in base::colnames(seurat_data@meta.data)){
    if (base::is.factor(seurat_data@meta.data[[i]])){
      base::print(base::paste("Re-evaluating levels for a factor column", i))
      base::print(
        base::paste(
          "before:", base::paste(base::levels(seurat_data@meta.data[[i]]), collapse=", ")
        )
      )
      seurat_data@meta.data[[i]] <- base::droplevels(seurat_data@meta.data[[i]])  # need to drop levels of the removed values
      base::print(
        base::paste(
          "after:", base::paste(base::levels(seurat_data@meta.data[[i]]), collapse=", ")
        )
      )
    }
  }
  return (seurat_data)
}

opt = NULL

options(warn = -1)

merge_and_recall_peaks = function(option_list)
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
  
  #### Read filtered object by doublets -----
  
  
  adata<-readRDS(file=opt$db_genotyped)
  
  # cat("adata_0\n")
  # cat(str(adata))
  # cat("\n")
  cat(sprintf(as.character(names(summary(as.factor(adata@meta.data$Assigned_GFPbc_integral))))))
  cat("\n")
  cat(sprintf(as.character(summary(as.factor(adata@meta.data$Assigned_GFPbc_integral)))))
  cat("\n")
  
  # Create a genotype factor ---------------------------------------------------------------------------------------
  
  

  adata@meta.data$Assigned_GFPgenotype_integral<-NA
  
  adata@meta.data$Assigned_GFPgenotype_integral[which(adata@meta.data$Assigned_GFPbc_integral%in%c('chrGFP_WTA','chrGFP_WTB','chrGFP_WTC'))]<-'wt'
  adata@meta.data$Assigned_GFPgenotype_integral[which(adata@meta.data$Assigned_GFPbc_integral%in%c('chrGFP_rs1','chrGFP_rs2','chrGFP_rs3'))]<-'CHEK2 T/T'
  adata@meta.data$Assigned_GFPgenotype_integral[which(adata@meta.data$Assigned_GFPbc_integral%in%c('chrGFP_R882H1','chrGFP_R882H2','chrGFP_R882H3'))]<-'DNMT3A R882H'
  adata@meta.data$Assigned_GFPgenotype_integral[which(adata@meta.data$Assigned_GFPbc_integral%in%c('chrGFP_rs_R882H1','chrGFP_rs_R882H2','chrGFP_rs_R882H3'))]<-'Double mutants'
  
  cat("adata@meta.data$Assigned_GFPgenotype_integral_PRE\n")
  cat(sprintf(as.character(names(summary(as.factor(adata@meta.data$Assigned_GFPgenotype_integral))))))
  cat("\n")
  cat(sprintf(as.character(summary(as.factor(adata@meta.data$Assigned_GFPgenotype_integral)))))
  cat("\n")
  
  adata@meta.data$Assigned_GFPgenotype_integral<-factor(adata@meta.data$Assigned_GFPgenotype_integral,
                                                        levels=c('wt','CHEK2 T/T','DNMT3A R882H','Double mutants'),
                                                        ordered=T)
  
  cat("adata@meta.data$Assigned_GFPgenotype_integral_POST\n")
  cat(sprintf(as.character(names(summary(as.factor(adata@meta.data$Assigned_GFPgenotype_integral))))))
  cat("\n")
  cat(sprintf(as.character(summary(as.factor(adata@meta.data$Assigned_GFPgenotype_integral)))))
  cat("\n")
  
  # Subset my object for cells with genotype assigned
  
  cat("SUBSET GENOTYPED\n")
  
  adata_geno<-subset(adata, Assigned_GFPgenotype_integral != "NA")
  
  cat("DROPLEVELS with new function\n")
  
  adata_geno<-refine_metadata_levels(adata_geno)
  
  
  
  
  

  
  #### Cluster without integration and check if any cluster looks specially bad for QC metrics -------------
  
  cat("RunUMAP RNA\n")
  
  DefaultAssay(adata_geno) <- 'RNA'
  
  adata_geno <- SCTransform(adata_geno, verbose = FALSE) 
  adata_geno <- RunPCA(adata_geno) 
  adata_geno <- RunUMAP(adata_geno, dims=1:50, reduction.name='umap.rna', reduction.key='rnaUMAP_')
  
  #### ATAC modality -------------
  # We exclude the first dimension as this is typically correlated with sequencing depth
  
  cat("RunUMAP ATAC\n")
  
  DefaultAssay(adata_geno) <- 'ATAC'

  adata_geno <- RunTFIDF(adata_geno)
  adata_geno <- FindTopFeatures(adata_geno, min.cutoff='q0')
  adata_geno <- RunSVD(adata_geno)
  adata_geno <- RunUMAP(adata_geno, reduction='lsi', dims=2:50, reduction.name='umap.atac', reduction.key='atacUMAP_')
  
  
  #### WNN ATAC+RNA modality -------------
  
  cat("RunUMAP wnn\n")
  
  adata_geno <- FindMultiModalNeighbors(adata_geno, reduction.list=list('pca', 'lsi'), dims.list=list(1:50, 2:50))
  adata_geno <- RunUMAP(adata_geno, nn.name='weighted.nn', reduction.name='umap.wnn', reduction.key='wnnUMAP_')
  adata_geno <- FindClusters(adata_geno, graph.name='wsnn', algorithm=4, resolution = 0.2, verbose=FALSE, method = "igraph")
  
  
  ###### SAVE -----
  
  setwd(out)
  
  saveRDS(adata_geno, file = 'merged_clusters_after_genotyping.rds')
  
  
  
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
    make_option(c("--db_genotyped"), type="character", default=NULL, 
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
  
  merge_and_recall_peaks(opt)
 

}


###########################################################################

system.time( main() )
  