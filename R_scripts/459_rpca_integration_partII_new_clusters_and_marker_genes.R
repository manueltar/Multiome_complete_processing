
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

rpca_integration = function(option_list)
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
  
  
  adata<-readRDS(file=opt$db_filt_clustered_QCed_cell_annotated_rpca_integrate)
  
  cat("adata_0\n")
  # cat(str(adata))
  cat("\n")
  
  #### Name the new graph -----
  
  cat("FindNeighbors\n")
  
  adata <- FindNeighbors(adata)
  
  #### Find new clusters at resolution 0.5 ----
  
  cat("FindClusters\n")
  
  
  adata <- FindClusters(adata, graph.name='integrated_nn', algorithm=4, resolution = 0.5, verbose=FALSE, method = "igraph")
  
  #### Prepare SCT for Find Markers see https://satijalab.org/seurat/reference/prepsctfindmarkers ----
  
  cat("PrepSCTFindMarkers\n")
  
  
  adata<-PrepSCTFindMarkers(adata, assay = "SCT", verbose = TRUE)
  
  #### Find marker genes ----
  
  cat("FindAllMarkers\n")
  
  
  adata.markers <- FindAllMarkers(adata, logfc.threshold = 0.25, test.use = "roc")
  
  cat("adata.markers_0\n")
  cat(str(adata.markers))
  cat("\n")
  
  

  ##### Saving ---------------
  
  cat("Saving\n")
  
  setwd(out)
  
  saveRDS(adata, file="merged_unprocessed_db_db_filt_clustered_QCed_cell_annotated_rpca_integrate_rpca_integrate_clustered.rds")
  
  saveRDS(adata.markers, file="merged_unprocessed_db_db_filt_clustered_QCed_cell_annotated_rpca_integrate_rpca_integrate_clustered_MARKER_GENES.rds")
  
  
  
  
 
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
    make_option(c("--db_filt_clustered_QCed_cell_annotated_rpca_integrate"), type="character", default=NULL, 
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
  
  rpca_integration(opt)
 

}


###########################################################################

system.time( main() )
  