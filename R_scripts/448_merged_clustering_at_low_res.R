
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
suppressMessages(library(EnsDb.Hsapiens.v86))
suppressMessages(library(dplyr))
suppressMessages(library(ggplot2))
suppressMessages(library(Matrix))
suppressMessages(library(data.table))
suppressMessages(library(ggpubr))
suppressMessages(library(ggplot2))
suppressMessages(library(scDblFinder))
suppressMessages(library("tidyr"))
suppressMessages(library("tibble"))
suppressMessages(library("biovizBase"))
suppressMessages(library("patchwork"))
suppressMessages(library(glmGamPoi))




opt = NULL

options(warn = -1)

cluster_at_low_res_and_exporth5ad = function(option_list)
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
  
  #### READ and transform res_param ----
  
  res_param = opt$res_param
  
  cat("res_param_\n")
  cat(sprintf(as.character(res_param)))
  cat("\n")
  
  #### Read filtered object by doublets -----
  
  
  adata2<-readRDS(file=opt$db_filt_clustered_QCed)
  
  cat("adata2_0\n")
  cat(str(adata2))
  cat("\n")
  
  
  
  #### Cluster without integration and check if any cluster looks specially bad for QC metrics -------------
  
  DefaultAssay(adata2) <- 'RNA'
  
  adata2 <- SCTransform(adata2, verbose = FALSE) 
  adata2 <- RunPCA(adata2) 
  adata2 <- RunUMAP(adata2, dims=1:50, reduction.name='umap.rna', reduction.key='rnaUMAP_')
  
  #### ATAC modality -------------
  # We exclude the first dimension as this is typically correlated with sequencing depth
  
  
  DefaultAssay(adata2) <- 'ATAC'

  adata2 <- RunTFIDF(adata2)
  adata2 <- FindTopFeatures(adata2, min.cutoff='q0')
  adata2 <- RunSVD(adata2)
  adata2 <- RunUMAP(adata2, reduction='lsi', dims=2:50, reduction.name='umap.atac', reduction.key='atacUMAP_')
  
  
  #### WNN ATAC+RNA modality -------------
  
  adata2 <- FindMultiModalNeighbors(adata2, reduction.list=list('pca', 'lsi'), dims.list=list(1:50, 2:50))
  adata2 <- RunUMAP(adata2, nn.name='weighted.nn', reduction.name='umap.wnn', reduction.key='wnnUMAP_')
  adata2 <- FindClusters(adata2, graph.name='wsnn', algorithm=4, resolution = res_param, verbose=FALSE, method = "igraph")
  
  
  ###### SAVE -----
  
  setwd(out)
  
  saveRDS(adata2, file = 'merged_unprocessed_db_filt_clustered_QCed_reclustered.rds')
  
  
  
  ##### Reduce the Seurat object to h5ad with RNA counts corrected by CellBender It doesn't work with CellBender corrected counts------
  
  DefaultAssay(adata2)<-'RNA'
  RNA_only<-DietSeurat(adata2, assays = "RNA")
  
  ###### SAVE -----
  
  setwd(out)
  
  SaveH5Seurat(RNA_only, filename = "merged_unprocessed_db_filt_clustered_QCed_reclustered.h5Seurat")
  Convert("merged_unprocessed_db_filt_clustered_QCed_reclustered.h5Seurat", dest = "h5ad")
  
  
  
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
    make_option(c("--db_filt_clustered_QCed"), type="character", default=NULL, 
                metavar="type", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--res_param"), type="numeric", default=NULL, 
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
  
  cluster_at_low_res_and_exporth5ad(opt)
 

}


###########################################################################

system.time( main() )
  