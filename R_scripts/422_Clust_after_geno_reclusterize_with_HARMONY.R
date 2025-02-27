

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
suppressMessages(library(harmony))



opt = NULL

options(warn = -1)



Harmony_recluster = function(option_list)
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
  
 
  
  ##### Read Seurat object----
  
  adata<-readRDS(file=opt$merged_clusters_after_genotyping)
  
  
  # RNA analysis
  
  cat("------------------------------------------------------------>Harmony_RNA\n")
  DefaultAssay(adata) <- 'RNA'
  adata <- SCTransform(adata, verbose=FALSE)
  adata <- RunPCA(adata)
  adata <- RunHarmony(adata, group.by.vars='orig.ident', assay.use='SCT', reduction.save='harmony.rna')
  adata <- RunUMAP(adata, dims=1:50, reduction='harmony.rna', reduction.name='umap.rna', reduction.key='rnaUMAP_')

  cat("------------------------------------------------------------>Harmony_ATAC\n")
  
  # ATAC analysis
  # We exclude the first dimension as this is typically correlated with sequencing depth
  DefaultAssay(adata) <- 'ATAC'
  adata <- RunTFIDF(adata)
  adata <- FindTopFeatures(adata, min.cutoff='q0', verbose=FALSE)
  adata <- RunSVD(adata)
  cat("------------------------------------------------------------>Harmony_ATAC_1\n")
  
  hm_atac <- HarmonyMatrix(Embeddings(adata, reduction='lsi'), adata$orig.ident, do_pca=FALSE)
  cat("------------------------------------------------------------>Harmony_ATAC_2\n")
  
  adata[['harmony.atac']] <- CreateDimReducObject(embeddings=hm_atac, key='atac_', assay='ATAC')
  
  cat("------------------------------------------------------------>Harmony_ATAC_3\n")
  
  adata <- RunUMAP(adata, dims=2:50, reduction='harmony.atac', reduction.name='umap.atac', reduction.key='atacUMAP_')

  cat("------------------------------------------------------------>Harmony_WNN\n")
  
  
  adata <- FindMultiModalNeighbors(adata, reduction.list=list('harmony.rna', 'harmony.atac'), dims.list=list(1:50, 2:50))
  adata <- RunUMAP(adata, nn.name='weighted.nn', reduction.name='umap.wnn', reduction.key='wnnUMAP_')
  adata <- FindClusters(adata, graph.name='wsnn', algorithm=4,  resolution = 0.2, verbose=FALSE)
  
  
  
  
 
 
  ##### SAVE RDS -------------------------
  
  
  setwd(out)
  
  
  saveRDS(adata, file="merged_clusters_after_genotyping_harmony_integrated.rds")
  
  
  
  
  
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
    make_option(c("--merged_clusters_after_genotyping"), type="character", default=NULL, 
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
  
  Harmony_recluster(opt)
 

}


###########################################################################

system.time( main() )
  