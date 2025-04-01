
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

create_the_premerged_Seurat_object = function(option_list)
{
  
  opt_in = option_list
  opt <<- option_list
  
  cat("All options:\n")
  printList(opt)
  
  
  #### READ and transform sample_name ----
  
  sample_name = opt$sample_name
  
  cat("sample_name_\n")
  cat(sprintf(as.character(sample_name)))
  cat("\n")
  
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
  
  path_processing_outputs = opt$path_processing_outputs
  
  cat("path_processing_outputs_0\n")
  cat(sprintf(as.character(path_processing_outputs)))
  cat("\n")
  
  intermediate_dir = opt$intermediate_dir
  
  cat("intermediate_dir_0\n")
  cat(sprintf(as.character(intermediate_dir)))
  cat("\n")
  
  snATAC_dir = opt$snATAC_dir
  
  cat("snATAC_dir_0\n")
  cat(sprintf(as.character(snATAC_dir)))
  cat("\n")
  
  crange_dir = opt$crange_dir
  
  cat("crange_dir_0\n")
  cat(sprintf(as.character(crange_dir)))
  cat("\n")
  
  premerge_dir = opt$premerge_dir
  
  cat("premerge_dir_0\n")
  cat(sprintf(as.character(premerge_dir)))
  cat("\n")
  
  fragfile = opt$fragfile
  
  cat("fragfile_0\n")
  cat(sprintf(as.character(fragfile)))
  cat("\n")
  
  ############### READ 5 kB ATAC windows per sample ------------------------
  
  filename<-paste(sample_name, '_snATAC_pipeline_job.long_fmt_mtx.txt.gz', sep='')
  
  cat("filename_\n")
  cat(sprintf(as.character(filename)))
  cat("\n")
  
  atac_lfmtx = read.table(paste0(snATAC_dir, filename))
  
  cat("atac_lfmtx_0\n")
  cat(str(atac_lfmtx))
  cat("\n")
  
  #### load preliminary_filtered seurat object ------------------------
  
  
  adata <- readRDS(file = opt$preliminary_filtered)
  
  # cat("adata_0\n")
  # cat(str(adata))
  # cat("\n")
  
  
  #### Read in cell bender counts -- filter them for previous bcs -- uses it as the main RNA modality ------------------
  
  cb = Read10X_h5(file.path(path_processing_outputs, 'cellbender_gex_seurat.h5'))
  cb_counts   <- cb$'Gene Expression'
  cb_counts   <- cb_counts[,colnames(adata)]
  
  cat("cb_counts_0\n")
  cat(str(cb_counts))
  cat("\n")
  
  adata2 = CreateSeuratObject(counts = cb_counts)
  
  adata2@meta.data$orig.ident = sample_name
  
  adata2[['percent.mt']] <- PercentageFeatureSet(adata2, pattern = '^MT-')
  
  cat("adata2_0\n")
  # cat(str(adata2))
  cat("\n")
  
  #add in previous raw RNA data as another assay (RNA_raw) for comparison
  
  DefaultAssay(adata) <- 'RNA'
  raw_rna <-  GetAssayData(object = adata, slot = "counts")
  raw_rna_assay <- CreateAssayObject(counts = raw_rna)
  adata2[['RNA_raw']] <- raw_rna_assay
  
  
  cat("adata2_1\n")
  # cat(str(adata2))
  cat("\n")
  
  
  
  ####### Add ATAC modality using previously generated 5kb windows matrix ---------------------------------
  
  FALSE %in% (colnames(adata2) %in%  atac_lfmtx$V1)
  FALSE %in% ( atac_lfmtx$V1 %in% colnames(adata2))
  
  atac_lfmtx$V1 <- factor(atac_lfmtx$V1, levels=colnames(adata2))
  reordered_lfm <- atac_lfmtx[order(atac_lfmtx$V1),]
  
  cat("atac_lfmtx_1\n")
  cat(str(atac_lfmtx))
  cat("\n")
  
  atac_sm <- with(reordered_lfm,
                  sparseMatrix(i=as.numeric(as.factor(V2)), j=as.numeric(V1),
                               x=V3, dimnames=list(levels(as.factor(V2)), levels(V1))))
  
  
  cat("atac_sm_0\n")
  cat(str(atac_sm))
  cat("\n")
  
  
  #create the new chromatin assay object and add to Seurat object
  
  atac_sm       <- atac_sm[,colnames(adata2)]
  grange.counts <- StringToGRanges(rownames(atac_sm), sep = c(':', '-'))
  grange.use    <- seqnames(grange.counts) %in% standardChromosomes(grange.counts)
  atac_sm       <- atac_sm[as.vector(grange.use), ]
  suppressMessages(annotations <- GetGRangesFromEnsDb(ensdb=EnsDb.Hsapiens.v86))
  seqlevelsStyle(annotations)  <- 'UCSC'
  genome(annotations)          <- 'hg38'
  
  cat("atac_sm_1\n")
  cat(str(atac_sm))
  cat("\n")
  
  suppressWarnings(chrom_assay <- CreateChromatinAssay(counts=atac_sm, sep=c(':', '-'),
                                                       genome='hg38', fragments=fragfile,
                                                       min.cells=-1, min.features=-1,
                                                       annotation=annotations))
  adata2[['ATAC']] <- chrom_assay
  
  invisible(gc())
  
  #### Compute new  metrics for cellbender rna and 5 kb windows atac -----------------------
  
  qc <- read.table(file.path(crange_dir, 'per_barcode_metrics.csv'), sep=',', header=TRUE, stringsAsFactors=1)
  qc <- as.data.frame(qc)
  rownames(qc) <- qc$gex_barcode
  qc <- qc[Cells(adata2), 6:length(colnames(qc))]
  adata2 <- AddMetaData(adata2, qc)
  
  DefaultAssay(adata2) <- 'ATAC'
  adata2 <- TSSEnrichment(adata2)
  
  
  ############ add previous metadata ------------------------------------------
  
  old_meta = adata@meta.data
  
  colkeep = c('scDblFinder.class','scDblFinder.score','scDblFinder.weighted','scDblFinder.cxds_score',
              'scDblFinder.class_atac','scDblFinder.score_atac','scDblFinder.weighted_atac','scDblFinder.cxds_score_atac',
              'No_assigned_GFPbc','Assigned_GFPbc','Assigned_GFPgenotype','DBL_comb')
  
  adata2@meta.data = cbind(adata2@meta.data,old_meta[,colkeep])
  
  ################# Read in Amulet barcodes and add to metadata ------------------------------
  
  amures           = read.table(file.path(intermediate_dir, "Amulet_selected_bc.tsv"))
  colnames(amures) = paste0("amulet_",colnames(amures))
  amures           = amures[rownames(adata2@meta.data), ]
  adata2@meta.data = cbind (adata2@meta.data, amures)
  adata2@meta.data$doublet_amulet = adata2@meta.data$amulet_q.value <0.05
  
  
  metadata_check<-adata2[[]]
  
  cat("metadata_check_0\n")
  cat(str(metadata_check))
  cat("\n")
  
  ################ CLUSTERIZATION ----------------------------------------------------------------------------
  
  adata2 <- adata2[, unname(which( colSums(GetAssayData(adata2, slot = "counts", assay = "RNA"))!=0))]
  
  ### Analyze and cluster
  # RNA analysis
  DefaultAssay(adata2) <- 'RNA'
  adata2 <- SCTransform(adata2, verbose = FALSE) 
  adata2 <- RunPCA(adata2) 
  adata2 <- RunUMAP(adata2, dims=1:50, reduction.name='umap.rna', reduction.key='rnaUMAP_')
  
  
  
  
  
  # ATAC analysis
  # We exclude the first dimension as this is typically correlated with sequencing depth
  DefaultAssay(adata2) <- 'ATAC'
  adata2 <- RunTFIDF(adata2)
  adata2 <- FindTopFeatures(adata2, min.cutoff='q0')
  adata2 <- RunSVD(adata2)
  adata2 <- RunUMAP(adata2, reduction='lsi', dims=2:50, reduction.name='umap.atac', reduction.key='atacUMAP_')
  
  # Multimodal analysis
  adata2 <- FindMultiModalNeighbors(adata2, reduction.list=list('pca', 'lsi'), dims.list=list(1:50, 2:50))
  adata2 <- RunUMAP(adata2, nn.name='weighted.nn', reduction.name='umap.wnn', reduction.key='wnnUMAP_')
  adata2 <- FindClusters(adata2, graph.name='wsnn', algorithm=4, resolution = .2, verbose=FALSE)
  
  
  
  
  
  
  
  ###### SAVE -----
  
  saveRDS(adata2, file = file.path(premerge_dir,'pre_merged.rds'))
  
  
  
}

graph_function = function(option_list)
{
  
  
  #### READ and transform sample_name ----
  
  sample_name = opt$sample_name
  
  cat("sample_name_\n")
  cat(sprintf(as.character(sample_name)))
  cat("\n")
  
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
  
  path_processing_outputs = opt$path_processing_outputs
  
  cat("path_processing_outputs_0\n")
  cat(sprintf(as.character(path_processing_outputs)))
  cat("\n")
  
  intermediate_dir = opt$intermediate_dir
  
  cat("intermediate_dir_0\n")
  cat(sprintf(as.character(intermediate_dir)))
  cat("\n")
  
  snATAC_dir = opt$snATAC_dir
  
  cat("snATAC_dir_0\n")
  cat(sprintf(as.character(snATAC_dir)))
  cat("\n")
  
  crange_dir = opt$crange_dir
  
  cat("crange_dir_0\n")
  cat(sprintf(as.character(crange_dir)))
  cat("\n")
  
  premerge_dir = opt$premerge_dir
  
  cat("premerge_dir_0\n")
  cat(sprintf(as.character(premerge_dir)))
  cat("\n")
  
  ##### Read the preprocessed object ------
  
  setwd(premerge_dir)
  
  
  
  adata2<-readRDS(file="pre_merged.rds")
  
  # cat("adata_0\n")
  # cat(str(adata2))
  # cat("\n")
  
  metadata_adata2<-adata2[[]]
  
  cat("metadata_adata2_0\n")
  cat(str(metadata_adata2))
  cat("\n")
  
  #### adata_assigned -------
  
  adata_assigned = subset(adata2, subset = No_assigned_GFPbc ==1)
  
  
  ########### Plot Intermediate_metrics ---------------
  
  cat("Plot Intermediate_metrics\n")
  
  sample_color = as.numeric(substring(sample_name, 7))
  ncells = read.csv(file.path(crange_dir, "summary.csv"))$Estimated.number.of.cells
  
  
  
  setwd(premerge_dir)
  
  Idents(adata2) <- "orig.ident"
  
  # options(repr.plot.width = 10, repr.plot.height = 5)
  # png(file.path(premerge_dir,'Intermediate_metrics.png'), width =800, height = 400)
  #   VlnPlot(adata2, features = c("nCount_ATAC", "nCount_RNA", "percent.mt",'TSS.enrichment'),  
  #           ncol = 4, cols=sample_color,
  #           log = TRUE, pt.size = 0)+ NoLegend()
  # dev.off()
  
  
  
  
  Vln_file<-VlnPlot(adata2, features = c("nCount_ATAC", "nCount_RNA", "percent.mt",'TSS.enrichment'),  
                    ncol = 4, cols=sample_color,
                    log = TRUE, pt.size = 0)+ NoLegend()
  
  saveRDS(Vln_file,file="Intermediate_metrics_graph.rds")
  
  
  # graph <- Vln_file +
  #   plot_layout( guides = "collect") & theme(text = element_text(size = 8, family = "Arial"))
  # 
  # # ggsave(p, filename = "plot.pdf", device=cairo_pdf, height=20, width=16, unit="cm", path=NULL, scale=1)
  # # ggsave(p, filename = "plot.svg", device=svg, height=20, width=16, unit="cm", path=NULL, scale=1)
  # 
  # ggsave(graph, filename = paste("Intermediate_metrics","_plot.pdf", sep=''), device=cairo_pdf, height=20, width=16, unit="cm", path=NULL, scale=1)
  
  
  # invisible(gc())
  
  ########### Plots of UMAPS and stats ---------------
  
  cat("Plots of UMAPS and stats\n")
  
  
  
  setwd(premerge_dir)
  
  options(repr.plot.width = 12, repr.plot.height = 5)
  png(file.path(premerge_dir,'Intermediate_UMAPs_clusters.png'), width =1000, height = 400)
  p1 <- DimPlot(adata2, reduction = "umap.rna", group.by = "seurat_clusters", label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle("RNA")
  p2 <- DimPlot(adata2, reduction = "umap.atac", group.by = "seurat_clusters", label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle("ATAC")
  p3 <- DimPlot(adata2, reduction = "umap.wnn", group.by = "seurat_clusters", label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle("WNN")
  p1 + p2 + p3 & NoLegend() & theme(plot.title = element_text(hjust = 0.5))
  dev.off()
  
  graph<-p1 + p2 + p3 & NoLegend() & theme(plot.title = element_text(hjust = 0.5))
  
  saveRDS(graph,file="UMAPS_and_stats_graph.rds")
  # 
  
  
  
  ########### Plot Intermediate_metrics and marker genes ---------------
  
  cat("Plot Intermediate_metrics and marker genes\n")
  
  
  DefaultAssay(adata2) <- 'SCT'
  
  
  
  ## Graph 
  
  DefaultAssay(adata) <- 'SCT'
  
  
  p5 <- FeaturePlot(adata, features = c('SOX2','POU5F1','EPCAM'),
                    reduction = 'umap.wnn', 
                    cols = c("lightgrey","darkgreen"), ncol = 3)
  
  p5_B <- FeaturePlot(adata, features = c("KDR",'CDH5','CD34'),
                      reduction = 'umap.wnn', 
                      cols = c("lightgrey","darkgreen"), ncol = 3)
  p5_1 <- FeaturePlot(adata, features = c('GYPA','HBA2','HBZ'),
                      reduction = 'umap.wnn', 
                      cols = c("lightgrey","darkgreen"), ncol = 3)
  
  p5_3 <- FeaturePlot(adata, features = c('ITGA2B','GP1BA', 'GP1BB'),
                      reduction = 'umap.wnn', 
                      cols = c("lightgrey","darkgreen"), ncol = 3)
  
  p_6 <- FeaturePlot(adata, features = c('TUBB1','CDKN2D','MAPK14'),
                     reduction = 'umap.wnn', 
                     cols = c("lightgrey","darkgreen"), ncol = 3)
  
  
  
  setwd(premerge_dir)
  
  options(repr.plot.width = 14, repr.plot.height = 8)
  png(file.path(premerge_dir,'Intermediate_UMAPs_qc.png'), width =1500, height = 1800)
  p5 / p5_B / p5_1 / p5_3 / p_6
  dev.off()
  
  
  graph<-p5 / p5_B / p5_1 / p5_3 / p_6
  
  saveRDS(graph,file="marker_genes_graph.rds")
  
  
  ### Plot Intermediate_UMAPs_doublets.png
  
  cat("Plot Intermediate_UMAPs_doublets.png\n")
  
  
  setwd(premerge_dir)
  
  png(file.path(premerge_dir,'Intermediate_UMAPs_doublets.png'), width =1000, height = 300)
  
  p1<- DimPlot(adata2, group.by = c("doublet_amulet"), reduction = 'umap.atac', )+ scale_color_manual(values = c("gray", "blue"))
  p2<- DimPlot(adata2, group.by = c("scDblFinder.class_atac"), reduction = 'umap.atac', )+ scale_color_manual(values = c("gray", "red4"))
  p3<- DimPlot(adata2, group.by = c("scDblFinder.class"), reduction = 'umap.atac', )+ scale_color_manual(values = c("gray", "green4"))
  
  
  p1 + p2 + p3
  dev.off()
  
  graph<-p1 + p2 +p3
  
  
  saveRDS(graph,file="Premerge_doublets.rds")
  
  #### Plot Intermediate_UMAPs_unique_GFP_bc
  
  cat("Plot Intermediate_UMAPs_unique_GFP_bc\n")
  
  
  setwd(premerge_dir)
  
  png(file.path(premerge_dir,'Intermediate_UMAPs_unique_GFP_bc.png'), width =1000, height = 600)
  p1<- DimPlot(adata_assigned, reduction = "umap.rna", group.by = "Assigned_GFPbc", label.size = 2.5, repel = TRUE) + ggtitle("RNA") 
  p2<- DimPlot(adata_assigned, reduction = "umap.atac", group.by = "Assigned_GFPbc", label.size = 2.5, repel = TRUE) + ggtitle("ATAC")
  p3<- DimPlot(adata_assigned, reduction = "umap.wnn", group.by = "Assigned_GFPbc", label.size = 2.5, repel = TRUE) + ggtitle("WNN")
  p4<- DimPlot(adata_assigned, reduction = "umap.wnn", group.by = "Assigned_GFPgenotype", label.size = 2.5, repel = TRUE) + ggtitle("WNN")
  
  p1 + p2 +p3 +p4 
  dev.off()
  
  graph<-p1 + p2 +p3 +p4 
  
  
  saveRDS(graph,file="Intermediate_UMAPs_unique_GFP_bc_graph.rds")
  
}

patchwork_function = function(option_list)
{
  opt_in = option_list
  opt <<- option_list
  
  cat("All options:\n")
  printList(opt)
  
  
  #### READ and transform sample_name ----
  
  sample_name = opt$sample_name
  
  cat("sample_name_\n")
  cat(sprintf(as.character(sample_name)))
  cat("\n")
  
  
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
  
  path_processing_outputs = opt$path_processing_outputs
  
  cat("path_processing_outputs_0\n")
  cat(sprintf(as.character(path_processing_outputs)))
  cat("\n")
  
  intermediate_dir = opt$intermediate_dir
  
  cat("intermediate_dir_0\n")
  cat(sprintf(as.character(intermediate_dir)))
  cat("\n")
  
  snATAC_dir = opt$snATAC_dir
  
  cat("snATAC_dir_0\n")
  cat(sprintf(as.character(snATAC_dir)))
  cat("\n")
  
  crange_dir = opt$crange_dir
  
  cat("crange_dir_0\n")
  cat(sprintf(as.character(crange_dir)))
  cat("\n")
  
  premerge_dir = opt$premerge_dir
  
  cat("premerge_dir_0\n")
  cat(sprintf(as.character(premerge_dir)))
  cat("\n")
  
  ##### Read the preprocessed object ------
  
  setwd(premerge_dir)
  
  
  
  ############## Retrieve graphs  -----------------
  
  
  file_list <- list.files(path=premerge_dir, include.dirs = FALSE)
  
  
  cat("file_list\n")
  cat(str(file_list))
  cat("\n")
  
  
  indexes_sel <- grep("_graph\\.rds",file_list)
  
  file_list_sel <- as.data.frame(file_list[indexes_sel], stringsAsFactors=F)
  colnames(file_list_sel)<-"file"
  
  cat("file_list_sel_0\n")
  cat(str(file_list_sel))
  cat("\n")
  
  file_list_sel$graph<-gsub("_graph\\.rds","",file_list_sel$file)
  
  
  cat("file_list_sel_1\n")
  cat(str(file_list_sel))
  cat("\n")  
  
  
  
  ############# LOOP --------------------------
  
  List_RESULTS<-list()
  
  
  # Results_DEF<-data.frame()
  
  DEBUG<-0
  
  for(i in 1:dim(file_list_sel)[1]){
    
    setwd(premerge_dir)
    
    read_file_sel<-file_list_sel$file[i]
    graph_sel<-file_list_sel$graph[i]
    
    
    cat("-------------------------------------------------------->\t")
    cat(sprintf(as.character(i)))
    cat("\t")
    cat(sprintf(as.character(read_file_sel)))
    cat("\t")
    cat(sprintf(as.character(graph_sel)))
    cat("\n")
    
    graph<-readRDS(file=read_file_sel)
    
    
    setwd(premerge_dir)
    
    p <- graph +
      plot_layout( guides = "collect") & theme(text = element_text(size = 8, family = "Arial"))
    
    # ggsave(p, filename = "plot.pdf", device=cairo_pdf, height=20, width=16, unit="cm", path=NULL, scale=1)
    # ggsave(p, filename = "plot.svg", device=svg, height=20, width=16, unit="cm", path=NULL, scale=1)
    
    ggsave(graph, filename = paste(graph_sel,"_plot.pdf", sep=''), device=cairo_pdf, height=40, width=32, unit="cm", path=NULL, scale=1)
    
    unlink(read_file_sel)
    
    List_RESULTS[[i]]<-graph
    
  }#i in 1:dim(file_list_sel)[1]
  
  
  # cat("List_RESULTS_0\n")
  # cat(str(List_RESULTS))
  # cat("\n")
  
  
  ########################## SAVE PLOTS -------------------------------------------------
  
  # setwd(premerge_dir)
  # 
  # p <- ((List_RESULTS[[1]] / List_RESULTS[[2]] / List_RESULTS[[3]]/ List_RESULTS[[4]]/ List_RESULTS[[5]])) +
  #   plot_annotation(tag_levels = "a") +
  #   plot_layout( guides = "collect") & theme(text = element_text(size = 8, family = "Arial"))
  # 
  # ggsave(p, filename = "plot.pdf", device=cairo_pdf, height=20, width=16, unit="cm", path=NULL, scale=1)
  # ggsave(p, filename = "plot.svg", device=svg, height=20, width=16, unit="cm", path=NULL, scale=1)
  
  
}

higher_resolution_clustering_in_relation_to_QC_metrics_graph_function = function(option_list)
{
  
  
  #### READ and transform sample_name ----
  
  sample_name = opt$sample_name
  
  cat("sample_name_\n")
  cat(sprintf(as.character(sample_name)))
  cat("\n")
  
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
  
  path_processing_outputs = opt$path_processing_outputs
  
  cat("path_processing_outputs_0\n")
  cat(sprintf(as.character(path_processing_outputs)))
  cat("\n")
  
  intermediate_dir = opt$intermediate_dir
  
  cat("intermediate_dir_0\n")
  cat(sprintf(as.character(intermediate_dir)))
  cat("\n")
  
  snATAC_dir = opt$snATAC_dir
  
  cat("snATAC_dir_0\n")
  cat(sprintf(as.character(snATAC_dir)))
  cat("\n")
  
  crange_dir = opt$crange_dir
  
  cat("crange_dir_0\n")
  cat(sprintf(as.character(crange_dir)))
  cat("\n")
  
  premerge_dir = opt$premerge_dir
  
  cat("premerge_dir_0\n")
  cat(sprintf(as.character(premerge_dir)))
  cat("\n")
  
  ##### Read the preprocessed object ------
  
  setwd(premerge_dir)
  
  
  
  adata2<-readRDS(file="pre_merged.rds")
  
  # cat("adata_0\n")
  # cat(str(adata2))
  # cat("\n")
  
  adata2 <- FindClusters(adata2, graph.name='wsnn', algorithm=4, resolution = 1, verbose=FALSE)
  
  p1 <- VlnPlot(adata2, features='nCount_SCT', group.by='seurat_clusters', log=TRUE) + geom_boxplot(width=.6, fill='white', alpha=.6) + geom_hline(yintercept=median(adata2$nCount_SCT), linetype='dashed')
  p2 <- VlnPlot(adata2, features='nFeature_SCT', group.by='seurat_clusters', log=TRUE) + geom_boxplot(width=.6, fill='white', alpha=.6) + geom_hline(yintercept=median(adata2$nFeature_SCT), linetype='dashed')
  p3 <- VlnPlot(adata2, features='nCount_ATAC', group.by='seurat_clusters', log=TRUE) + geom_boxplot(width=.6, fill='white', alpha=.6) + geom_hline(yintercept=median(adata2$nCount_ATAC), linetype='dashed')
  p4 <- VlnPlot(adata2, features='nFeature_ATAC', group.by='seurat_clusters', log=TRUE) + geom_boxplot(width=.6, fill='white', alpha=.6) + geom_hline(yintercept=median(adata2$nFeature_ATAC), linetype='dashed')
  p5 <- VlnPlot(adata2, features='nCount_RNA', group.by='seurat_clusters', log=TRUE) + geom_boxplot(width=.6, fill='white', alpha=.6) + geom_hline(yintercept=median(adata2$nCount_SCT), linetype='dashed')
  p6 <- VlnPlot(adata2, features='nFeature_RNA', group.by='seurat_clusters', log=TRUE) + geom_boxplot(width=.6, fill='white', alpha=.6) + geom_hline(yintercept=median(adata2$nFeature_SCT), linetype='dashed')
  
  
  #options(repr.plot.width=24, repr.plot.height=8)
  #ggarrange(p1, p3, p5, p2, p4, p6, ncol = 3, nrow = 2)
  
  
  p9 <- VlnPlot(adata2, features='TSS.enrichment', group.by='seurat_clusters', log=TRUE) + geom_boxplot(width=.6, fill='white', alpha=.6) + geom_hline(yintercept=median(adata2$TSS.enrichment), linetype='dashed')
  p10 <- VlnPlot(adata2, features='percent.mt', group.by='seurat_clusters', log=TRUE) + geom_boxplot(width=.6, fill='white', alpha=.6) + geom_hline(yintercept=median(adata2$percent.mt), linetype='dashed')
  p11 <- VlnPlot(adata2, features='amulet_nFrags', group.by='seurat_clusters', log=TRUE) + geom_boxplot(width=.6, fill='white', alpha=.6) + geom_hline(yintercept=median(adata2$amulet_nFrags), linetype='dashed')
  p12 <- VlnPlot(adata2, features='scDblFinder.score', group.by='seurat_clusters') + geom_boxplot(width=.6, fill='white', alpha=.6) + geom_hline(yintercept=median(adata2$scDblFinder.score), linetype='dashed')
  p13 <- VlnPlot(adata2, features='scDblFinder.score_atac', group.by='seurat_clusters') + geom_boxplot(width=.6, fill='white', alpha=.6) + geom_hline(yintercept=median(adata2$scDblFinder.score_atac), linetype='dashed')
  p14 <- VlnPlot(adata2, features='amulet_q.value', group.by='seurat_clusters') + geom_boxplot(width=.6, fill='white', alpha=.6) + geom_hline(yintercept=median(adata2$scDblFinder.score_atac), linetype='dashed')
  
  
  #options(repr.plot.width=24, repr.plot.height=8)
  #ggarrange(p9, p10, p11, p12, p13, p14, ncol = 3, nrow = 2)
  
  print(
    png(file.path(premerge_dir,'Violin_plots_QC_byCluster.png'), width =1200, height = 1200))
  ggarrange(p1& NoLegend() , p3& NoLegend(), p5& NoLegend(), p2& NoLegend(), p4& NoLegend(), p6& NoLegend(),
            p9& NoLegend(), p10& NoLegend(), p11& NoLegend(), p12& NoLegend(), p13& NoLegend(), p14 & NoLegend(), ncol = 3, nrow = 4)
  dev.off()
  
  print(
    png(file.path(premerge_dir,'Premerge_UMAPs_clusters_res1.png'), width =1000, height = 350))
  p1 <- DimPlot(adata2, reduction = "umap.rna", group.by = "seurat_clusters", label = TRUE,  repel = TRUE) + ggtitle("RNA")
  p2 <- DimPlot(adata2, reduction = "umap.atac", group.by = "seurat_clusters", label = TRUE, repel = TRUE) + ggtitle("ATAC")
  p3 <- DimPlot(adata2, reduction = "umap.wnn", group.by = "seurat_clusters", label = TRUE,  repel = TRUE) + ggtitle("WNN")
  p1 + p2 + p3 & NoLegend() & theme(plot.title = element_text(hjust = 0.5))
  dev.off()
  
  ### Any sample-specific clusters?
  
  p1 <- VlnPlot(adata2, features='nCount_SCT', group.by='Assigned_GFPbc', log=TRUE) + geom_boxplot(width=.6, fill='white', alpha=.6) + geom_hline(yintercept=median(adata2$nCount_SCT), linetype='dashed')
  p2 <- VlnPlot(adata2, features='nFeature_SCT', group.by='Assigned_GFPbc', log=TRUE) + geom_boxplot(width=.6, fill='white', alpha=.6) + geom_hline(yintercept=median(adata2$nFeature_SCT), linetype='dashed')
  p3 <- VlnPlot(adata2, features='nCount_ATAC', group.by='Assigned_GFPbc', log=TRUE) + geom_boxplot(width=.6, fill='white', alpha=.6) + geom_hline(yintercept=median(adata2$nCount_ATAC), linetype='dashed')
  p4 <- VlnPlot(adata2, features='nFeature_ATAC', group.by='Assigned_GFPbc', log=TRUE) + geom_boxplot(width=.6, fill='white', alpha=.6) + geom_hline(yintercept=median(adata2$nFeature_ATAC), linetype='dashed')
  p5 <- VlnPlot(adata2, features='nCount_RNA', group.by='Assigned_GFPbc', log=TRUE) + geom_boxplot(width=.6, fill='white', alpha=.6) + geom_hline(yintercept=median(adata2$nCount_SCT), linetype='dashed')
  p6 <- VlnPlot(adata2, features='nFeature_RNA', group.by='Assigned_GFPbc', log=TRUE) + geom_boxplot(width=.6, fill='white', alpha=.6) + geom_hline(yintercept=median(adata2$nFeature_SCT), linetype='dashed')
  
  
  #options(repr.plot.width=24, repr.plot.height=8)
  #ggarrange(p1, p3, p5, p2, p4, p6, ncol = 3, nrow = 2)
  
  p9 <- VlnPlot(adata2, features='TSS.enrichment', group.by='Assigned_GFPbc', log=TRUE) + geom_boxplot(width=.6, fill='white', alpha=.6) + geom_hline(yintercept=median(adata2$TSS.enrichment), linetype='dashed')
  p10 <- VlnPlot(adata2, features='percent.mt', group.by='Assigned_GFPbc', log=TRUE) + geom_boxplot(width=.6, fill='white', alpha=.6) + geom_hline(yintercept=median(adata2$percent.mt), linetype='dashed')
  p11 <- VlnPlot(adata2, features='amulet_nFrags', group.by='Assigned_GFPbc', log=TRUE) + geom_boxplot(width=.6, fill='white', alpha=.6) + geom_hline(yintercept=median(adata2$amulet_nFrags), linetype='dashed')
  p12 <- VlnPlot(adata2, features='scDblFinder.score', group.by='Assigned_GFPbc') + geom_boxplot(width=.6, fill='white', alpha=.6) + geom_hline(yintercept=median(adata2$scDblFinder.score), linetype='dashed')
  p13 <- VlnPlot(adata2, features='scDblFinder.score_atac', group.by='Assigned_GFPbc') + geom_boxplot(width=.6, fill='white', alpha=.6) + geom_hline(yintercept=median(adata2$scDblFinder.score_atac), linetype='dashed')
  p14 <- VlnPlot(adata2, features='amulet_q.value', group.by='Assigned_GFPbc') + geom_boxplot(width=.6, fill='white', alpha=.6) + geom_hline(yintercept=median(adata2$scDblFinder.score_atac), linetype='dashed')
  
  
  #options(repr.plot.width=24, repr.plot.height=8)
  #ggarrange(p9, p10, p11, p12, p13, p14, ncol = 3, nrow = 2)
  
  
  qc_clus = data.frame(barcodes = rownames(adata2@meta.data), clusters = adata2@meta.data[,"seurat_clusters"])
  write.table(qc_clus, file.path(premerge_dir,'qc_clusters.tsv'), quote=F, row.names=FALSE)
  print(png(file.path(premerge_dir,'Violin_plots_QC_bysample.png'), width =1200, height = 1200))
  ggarrange(p1& NoLegend() , p3& NoLegend(), p5& NoLegend(), p2& NoLegend(), p4& NoLegend(), p6& NoLegend(),
            p9& NoLegend(), p10& NoLegend(), p11& NoLegend(), p12& NoLegend(), p13& NoLegend(), p14 & NoLegend(), ncol = 3, nrow = 4)
  dev.off()
  
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
    make_option(c("--sample_name"), type="character", default=NULL, 
                metavar="type", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--fragfile"), type="character", default=NULL, 
                metavar="type", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--path_processing_outputs"), type="character", default=NULL, 
                metavar="type", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--intermediate_dir"), type="character", default=NULL, 
                metavar="type", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--snATAC_dir"), type="character", default=NULL, 
                metavar="type", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--crange_dir"), type="character", default=NULL, 
                metavar="type", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--premerge_dir"), type="character", default=NULL, 
                metavar="type", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--preliminary_filtered"), type="character", default=NULL, 
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
  
  create_the_premerged_Seurat_object(opt)
  graph_function(opt)
  patchwork_function(opt)
  higher_resolution_clustering_in_relation_to_QC_metrics_graph_function(opt)
  
  
}


###########################################################################

system.time( main() )
