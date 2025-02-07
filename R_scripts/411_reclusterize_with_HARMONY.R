
.libPaths()
.libPaths(new = c("/home/manuel.tardaguila/conda_envs/multiome_QC/lib/R/library"))
.libPaths()
# sessionInfo()

Sys.setenv(RETICULATE_PYTHON="/home/manuel.tardaguila/conda_envs/multiome_QC/bin/python")
library(reticulate)
reticulate::use_python("/home/manuel.tardaguila/conda_envs/multiome_QC/bin/python")
reticulate::use_condaenv("/home/manuel.tardaguila/conda_envs/multiome_QC")
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
  
  path_graphs = paste(out,'graphs_FINAL_HARMONY','/',sep='')
  
  if (file.exists(path_graphs)){
    
    
  }else{
    
    dir.create(file.path(path_graphs))
    
  }#path_processing_outputs
  
  
  ##### Read Seurat object----
  
  adata3<-readRDS(file=opt$db_filt_clustered_minus_26_MACS2_peaks)
  
  
  # RNA analysis
  
  cat("------------------------------------------------------------>Harmony_RNA\n")
  DefaultAssay(adata3) <- 'RNA'
  adata3 <- SCTransform(adata3, verbose=FALSE)
  adata3 <- RunPCA(adata3)
  adata3 <- RunHarmony(adata3, group.by.vars='orig.ident', assay.use='SCT', reduction.save='harmony.rna')
  adata3 <- RunUMAP(adata3, dims=1:50, reduction='harmony.rna', reduction.name='umap.rna', reduction.key='rnaUMAP_')

  cat("------------------------------------------------------------>Harmony_ATAC\n")
  
  # ATAC analysis
  # We exclude the first dimension as this is typically correlated with sequencing depth
  DefaultAssay(adata3) <- 'ATAC'
  adata3 <- RunTFIDF(adata3)
  adata3 <- FindTopFeatures(adata3, min.cutoff='q0', verbose=FALSE)
  adata3 <- RunSVD(adata3)
  cat("------------------------------------------------------------>Harmony_ATAC_1\n")
  
  hm_atac <- HarmonyMatrix(Embeddings(adata3, reduction='lsi'), adata3$orig.ident, do_pca=FALSE)
  cat("------------------------------------------------------------>Harmony_ATAC_2\n")
  
  adata3[['harmony.atac']] <- CreateDimReducObject(embeddings=hm_atac, key='atac_', assay='ATAC')
  
  cat("------------------------------------------------------------>Harmony_ATAC_3\n")
  
  adata3 <- RunUMAP(adata3, dims=2:50, reduction='harmony.atac', reduction.name='umap.atac', reduction.key='atacUMAP_')

  cat("------------------------------------------------------------>Harmony_WNN\n")
  
  
  adata3 <- FindMultiModalNeighbors(adata3, reduction.list=list('harmony.rna', 'harmony.atac'), dims.list=list(1:50, 2:50))
  adata3 <- RunUMAP(adata3, nn.name='weighted.nn', reduction.name='umap.wnn', reduction.key='wnnUMAP_')
  adata3 <- FindClusters(adata3, graph.name='wsnn', algorithm=4,  resolution = 0.5, verbose=FALSE)
  
  
  
  
  ## Graph WNN by SampleID
  
  p3 <- DimPlot(adata3, reduction = "umap.wnn", group.by = "orig.ident", label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle("WNN")
  
  setwd(path_graphs)
  png(file='WNN_by_SampleID.png', width =500, height = 350)
  p3  & theme(plot.title = element_text(hjust = 0.5))
  dev.off()
  
  ## Graph WNN by Seurat cluster
  
  p3 <- DimPlot(adata3, reduction = "umap.wnn", group.by = "seurat_clusters", label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle("WNN")
  
  setwd(path_graphs)
  png(file='WNN_by_Seurat_cluster.png', width =500, height = 350)
  p3  & theme(plot.title = element_text(hjust = 0.5))
  dev.off()
  
  ## Graph WNN by GFPbc
  
  p3 <- DimPlot(adata3, reduction = "umap.wnn", group.by = "Assigned_GFPbc", label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle("WNN")
  
  setwd(path_graphs)
  png(file='WNN_by_Assigned_GFPbc.png', width =500, height = 350)
  p3  & theme(plot.title = element_text(hjust = 0.5))
  dev.off()
  
  ## Graph WNN by genotype
  
  p3 <- DimPlot(adata3, reduction = "umap.wnn", group.by = "Assigned_GFPgenotype", label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle("WNN")
  
  setwd(path_graphs)
  png(file='WNN_by_Assigned_GFPgenotype.png', width =500, height = 350)
  p3  & theme(plot.title = element_text(hjust = 0.5))
  dev.off()
  
  ## Graph WNN_Intermediate_QC_metrics
  
  p6 <- FeaturePlot(adata3, features = c("nCount_SCT", "nCount_RNA", "nCount_ATAC",'TSS.enrichment'), ncol = 4,
                    reduction = 'umap.wnn')
  
  
  
  setwd(path_graphs)
  png(file='WNN_Intermediate_QC_metrics.png', width =1000, height = 500)
  p6
  dev.off()
  
  p6 <- FeaturePlot(adata3, features = c("nFeature_SCT", "nFeature_RNA", "nFeature_ATAC",'percent.mt'), ncol = 4,
                    reduction = 'umap.wnn')
  
  
  
  setwd(path_graphs)
  png(file='WNN_Intermediate_QC_metrics_2.png', width =1000, height = 500)
  p6
  dev.off()
  
  ## Graph WNN by marker genes
  
  DefaultAssay(adata3) <- 'SCT'
  
  
  p5 <- FeaturePlot(adata3, features = c("NANOG",'SOX4', 'EOMES', 'TWIST1'),
                    reduction = 'umap.wnn', 
                    cols = c("lightgrey","darkgreen"), ncol = 4)
  p5_1 <- FeaturePlot(adata3, features = c('GYPA','HBA2','HBZ','ITGA2B'),
                      reduction = 'umap.wnn', 
                      cols = c("lightgrey","darkgreen"), ncol = 4)
  p5_2 <- FeaturePlot(adata3, features = c('GP1BA','KIF15','STIL','ANAPC7'),
                      reduction = 'umap.wnn', 
                      cols = c("lightgrey","darkgreen"), ncol = 4)
  
  setwd(path_graphs)
  png(file='WNN_marker_genes.png', width =1000, height = 500)
  p5 / p5_1 /p5_2
  dev.off()
  
  # Violin graphs to decide if a cluster is lowQuality
  
  p10 <- VlnPlot(adata3, features='percent.mt', group.by='seurat_clusters', log=TRUE) + geom_boxplot(width=.6, fill='white', alpha=.6) + geom_hline(yintercept=median(adata3$percent.mt), linetype='dashed')
  
  
  
  setwd(path_graphs)
  png(file='Vln_percent.mt.png', width =1000, height = 500)
  p10
  dev.off()
  
  ## nCount_SCT
  
  p10 <- VlnPlot(adata3, features='nCount_SCT', group.by='seurat_clusters', log=TRUE) + geom_boxplot(width=.6, fill='white', alpha=.6) + geom_hline(yintercept=median(adata3$nCount_SCT), linetype='dashed')
  
  
  
  setwd(path_graphs)
  png(file='Vln_nCount_SCT.png', width =1000, height = 500)
  p10
  dev.off()
  
  ## nFeature_SCT
  
  p10 <- VlnPlot(adata3, features='nFeature_SCT', group.by='seurat_clusters', log=TRUE) + geom_boxplot(width=.6, fill='white', alpha=.6) + geom_hline(yintercept=median(adata3$nFeature_SCT), linetype='dashed')
  
  
  
  setwd(path_graphs)
  png(file='Vln_nFeature_SCT.png', width =1000, height = 500)
  p10
  dev.off()
  
  ## nCount_ATAC
  
  p10 <- VlnPlot(adata3, features='nCount_ATAC', group.by='seurat_clusters', log=TRUE) + geom_boxplot(width=.6, fill='white', alpha=.6) + geom_hline(yintercept=median(adata3$nFeature_ATAC), linetype='dashed')
  
  
  setwd(path_graphs)
  png(file='Vln_nCount_ATAC.png', width =1000, height = 500)
  p10
  dev.off()
  
  
  
  ## nFeature_ATAC
  
  ####
  
  p10 <- VlnPlot(adata3, features='nFeature_ATAC', group.by='seurat_clusters', log=TRUE) + geom_boxplot(width=.6, fill='white', alpha=.6) + geom_hline(yintercept=median(adata3$nCount_RNA), linetype='dashed')
  
  setwd(path_graphs)
  png(file='Vln_nFeature_ATAC.png', width =1000, height = 500)
  p10
  dev.off()
  
  
  
  ## nCount_RNA
  
  ####
  
  p10 <- VlnPlot(adata3, features='nCount_RNA', group.by='seurat_clusters', log=TRUE) + geom_boxplot(width=.6, fill='white', alpha=.6) + geom_hline(yintercept=median(adata3$nCount_RNA), linetype='dashed')
  
  
  
  setwd(path_graphs)
  png(file='Vln_nCount_RNA.png', width =1000, height = 500)
  p10
  dev.off()
  
  
  ## nFeature_RNA
  
  
  ####
  
  p10 <- VlnPlot(adata3, features='nFeature_RNA', group.by='seurat_clusters', log=TRUE) + geom_boxplot(width=.6, fill='white', alpha=.6) + geom_hline(yintercept=median(adata3$nFeature_RNA), linetype='dashed')
  
  
  
  setwd(path_graphs)
  png(file='Vln_nFeature_RNA.png', width =1000, height = 500)
  p10
  dev.off()
  
  
  ## TSS.enrichment
  
  ####
  
  p10 <- VlnPlot(adata3, features='TSS.enrichment', group.by='seurat_clusters', log=TRUE) + geom_boxplot(width=.6, fill='white', alpha=.6) + geom_hline(yintercept=median(adata3$TSS.enrichment), linetype='dashed')
  
  
  
  setwd(path_graphs)
  png(file='Vln_TSS.enrichment.png', width =1000, height = 500)
  p10
  dev.off()
  
  ## amulet_nFrags
  
  ####
  
  p10 <- VlnPlot(adata3, features='amulet_nFrags', group.by='seurat_clusters', log=TRUE) + geom_boxplot(width=.6, fill='white', alpha=.6) + geom_hline(yintercept=median(adata3$amulet_nFrags), linetype='dashed')
  
  
  
  setwd(path_graphs)
  png(file='Vln_amulet_nFrags.png', width =1000, height = 500)
  p10
  dev.off()
  
  
  ## scDblFinder.score
  
  ####
  
  p10 <- VlnPlot(adata3, features='scDblFinder.score', group.by='seurat_clusters', log=TRUE) + geom_boxplot(width=.6, fill='white', alpha=.6) + geom_hline(yintercept=median(adata3$scDblFinder.score), linetype='dashed')
  
  
  
  setwd(path_graphs)
  png(file='Vln_scDblFinder.score.png', width =1000, height = 500)
  p10
  dev.off()
  
  
  ## scDblFinder.score_atac
  
  ####
  
  p10 <- VlnPlot(adata3, features='scDblFinder.score_atac', group.by='seurat_clusters', log=TRUE) + geom_boxplot(width=.6, fill='white', alpha=.6) + geom_hline(yintercept=median(adata3$scDblFinder.score_atac), linetype='dashed')
  
  
  
  setwd(path_graphs)
  png(file='Vln_scDblFinder.score_atac.png', width =1000, height = 500)
  p10
  dev.off()
  
 
  ##### SAVE RDS -------------------------
  
  
  setwd(out)
  
  
  saveRDS(adata3, file="merged_unprocessed_db_filt_clustered_minus_26_MACS2_peaks_HARMONY_clustered.rds")
  
  
  
  
  
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
    make_option(c("--db_filt_clustered_minus_26_MACS2_peaks"), type="character", default=NULL, 
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
  