
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




opt = NULL

options(warn = -1)

create_the_preliminary_Seurat_object = function(option_list)
{
  
  opt_in = option_list
  opt <<- option_list
  
  cat("All options:\n")
  printList(opt)
  
  
  Uniquely_genotyped_larry_barcodes_assignments<-readRDS(opt$Uniquely_genotyped_larry_barcodes_assignments)
  
  cat("Uniquely_genotyped_larry_barcodes_assignments_0\n")
  cat(str(Uniquely_genotyped_larry_barcodes_assignments))
  cat("\n")
  
  #### READ and transform Threshold_UMIS_per_cell ----
  
  Threshold_UMIS_per_cell = opt$Threshold_UMIS_per_cell
  
  cat("Threshold_UMIS_per_cell_\n")
  cat(sprintf(as.character(Threshold_UMIS_per_cell)))
  cat("\n")
  
  
  #### READ and transform sample_name ----
  
  sample_name = opt$sample_name
  
  cat("sample_name_\n")
  cat(sprintf(as.character(sample_name)))
  cat("\n")
  
  #### READ and transform master_path ----
  
  master_path = opt$master_path
  
  cat("master_path_\n")
  cat(sprintf(as.character(master_path)))
  cat("\n")
  
  sample_dir<-paste(master_path,sample_name,'/','outs','/',sep='')
  
  cat("sample_dir_\n")
  cat(sprintf(as.character(sample_dir)))
  cat("\n")
  
  #### READ and transform rna_min_features ----
  
  rna_min_features = opt$rna_min_features
  
  cat("rna_min_features_\n")
  cat(sprintf(as.character(rna_min_features)))
  cat("\n")
  
  #### READ and transform atac_min_fragments ----
  
  atac_min_fragments = opt$atac_min_fragments
  
  cat("atac_min_fragments_\n")
  cat(sprintf(as.character(atac_min_fragments)))
  cat("\n")
  
  #### READ and transform MITO_max ----
  
  MITO_max = opt$MITO_max
  
  cat("MITO_max_\n")
  cat(sprintf(as.character(MITO_max)))
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
  
  path_processing_outputs = paste(out,sample_name,'/',sep='')
  
  if (file.exists(path_processing_outputs)){
    
    
  }else{
    
    dir.create(file.path(path_processing_outputs))
    
  }#path_processing_outputs
  
  
  output_dir = paste(path_processing_outputs,'intermediate','/',sep='')
  
  if (file.exists(output_dir)){
    
    
  }else{
    
    dir.create(file.path(output_dir))
    
  }#output_dir
  
  
  
  
  ###### Read raw data -------
  
  inputdata.10x         <- Read10X_h5(file.path(sample_dir, 'raw_feature_bc_matrix.h5'))
  
  cat("inputdata.10x_0\n")
  cat(str(inputdata.10x))
  cat("\n")
  
  rna_counts            <- inputdata.10x$'Gene Expression'
  atac_counts           <- inputdata.10x$'Peaks'
  
  cat("atac_counts_0\n")
  cat(str(atac_counts))
  cat("\n")
  
  ## Create object with rna_counts
  
  adata                 <- CreateSeuratObject(counts=rna_counts)
  
  cat("adata_0\n")
  cat(str(adata))
  cat("\n")
  
  metadata_adata<-adata[[]]
  
  cat("metadata_adata_0\n")
  cat(str(metadata_adata))
  cat("\n")
  
  ### Add the percent.mt for genes that start with MT (mithochondrial genes)
  
  adata[['percent.mt']] <- PercentageFeatureSet(adata, pattern = '^MT-')
  
  
  ##### Create the filter for cells that have equal or more than rna_min_features genes -----
  
  stored_filters = c() # Define the object to store the CB that pass the subsquent filters
  
  stored_filters['total_bc'] = length(colnames(adata[["RNA"]])) # All the cells of the adata RNA layer
  
  
  adata_sub <- subset(x = adata,subset = nFeature_RNA >= as.numeric(rna_min_features)) # subset for cells with more or equal to 500 genes
  
  cat("adata_sub_0\n")
  cat(str(adata_sub))
  cat("\n")
  
  stored_filters['after_rna_minfeat'] = length(colnames(adata_sub[["RNA"]])) # Only the cells with more or equal to 500 genes
  
  cat("stored_filters_0\n")
  cat(str(stored_filters))
  cat("\n")
  
  ##### add cellranger QC metadata to adata_sub object --------------
  
  qc <- read.table(file.path(sample_dir, 'per_barcode_metrics.csv'), sep=',', header=TRUE, stringsAsFactors=1)
  qc           <- as.data.frame(qc)
  rownames(qc) <- qc$gex_barcode
  qc           <- qc[Cells(adata_sub), ]
  adata_sub    <- AddMetaData(adata_sub, qc)
  
  metadata_adata_sub<-adata_sub[[]]
  
  cat("metadata_adata_sub_0\n")
  cat(str(metadata_adata_sub))
  cat("\n")
  
  ####### Add in ATAC data for the barcodes that passed RNA filters ----------
  
  # Subset atac counts to cell barcodes that passed the rna_min_features filter
  
  atac_counts   <- atac_counts[,colnames(adata_sub)]
  
  cat("atac_counts_1\n")
  cat(str(atac_counts))
  cat("\n")
  
  # Make the format of the peaks compatible with frag.file
  
  grange.counts <- StringToGRanges(rownames(atac_counts), sep = c(':', '-'))
  grange.use    <- seqnames(grange.counts) %in% standardChromosomes(grange.counts)
  atac_counts   <- atac_counts[as.vector(grange.use), ]
  
  cat("atac_counts_2\n")
  cat(str(atac_counts))
  cat("\n")
  
  suppressMessages(annotations <- GetGRangesFromEnsDb(ensdb=EnsDb.Hsapiens.v86))
  seqlevelsStyle(annotations)  <- 'UCSC'
  genome(annotations)          <- 'hg38'
  
  
  ### Read frag.file
  
  frag.file <- file.path(sample_dir, 'atac_fragments.tsv.gz')
  
  cat("frag.file_0\n")
  cat(str(frag.file))
  cat("\n")
  
  ### Create the chrom_assay layer
  
  suppressWarnings(chrom_assay <- CreateChromatinAssay(counts=atac_counts, sep=c(':', '-'), 
                                                       genome='hg38', fragments=frag.file, 
                                                       min.cells=-1, min.features=-1, 
                                                       annotation=annotations))
  ### Add ATAC layer to adata_sub
  
  adata_sub[['ATAC']] <- chrom_assay
  
  invisible(gc()) ### ????
  
  cat("adata_sub_1\n")
  cat(str(adata_sub))
  cat("\n")
  
  ########## filter by min ATAC number fragments keep only cells that have 1000 or more ATAC peaks  ------------------
  
  adata_sub_atac <- subset( x = adata_sub, subset = atac_fragments >= as.numeric(atac_min_fragments))
  
  
  stored_filters['after_atac_minfrag'] = length(colnames(adata_sub_atac[["RNA"]]))
  
  cat("stored_filters_1\n")
  cat(str(stored_filters))
  cat("\n")
  
  ####### First doublet detection RNA with scDblFinder to be performed on lightly filtered data on the already filtered adata_sub_atac object ------
  
  DefaultAssay(adata_sub_atac) <- 'RNA'
  sce <- scDblFinder(GetAssayData(object = adata_sub_atac, slot = "counts"))
  sce_results = data.frame(SummarizedExperiment::colData(sce))
  adata_sub_atac@meta.data = cbind(adata_sub_atac@meta.data,sce_results)
  
  
  metadata_adata_sub_atac<-adata_sub_atac[[]]
  
  cat("metadata_adata_sub_atac_0\n")
  cat(str(metadata_adata_sub_atac))
  cat("\n")
  
  
  ### First doublet detection ATAC with scDblFinder to be performed on lightly filtered data on the already filtered adata_sub_atac object ----------------------------
  
  DefaultAssay(adata_sub_atac) <- 'ATAC'
  sce_atac <- scDblFinder(GetAssayData(object = adata_sub_atac, slot = "counts"), 
                          artificialDoublets=1, aggregateFeatures=TRUE, 
                          nfeatures=25, processing="normFeatures")
  
  
  sce_results_atac = data.frame(SummarizedExperiment::colData(sce_atac))
  colnames(sce_results_atac) = paste(colnames(sce_results_atac), "atac", sep="_")
  adata_sub_atac@meta.data = cbind(adata_sub_atac@meta.data,sce_results_atac)
  
  metadata_adata_sub_atac<-adata_sub_atac[[]]
  
  cat("metadata_adata_sub_atac_1\n")
  cat(str(metadata_adata_sub_atac))
  cat("\n")
  
  
  ### Add a metadata column that has both scDblFinder.class and scDblFinder.class_atac
  
  adata_sub_atac@meta.data$DBL_comb = paste("R",adata_sub_atac@meta.data$scDblFinder.class, 
                                   "A",adata_sub_atac@meta.data$scDblFinder.class_atac, sep=":") 
  
  metadata_adata_sub_atac<-adata_sub_atac[[]]
  
  
  cat("metadata_adata_sub_atac_2\n")
  cat(str(metadata_adata_sub_atac))
  cat("\n")
  
  ##### Remove multiplets and other cellranger excluded cells creating the adata_sub_multiplet from the adata_sub_atac ----------------------
  
  adata_sub_multiplet <- subset(
    x = adata_sub_atac,
    subset = excluded_reason != 1 
  )
  stored_filters['after_cr_multiplets'] = length(colnames(adata_sub_atac[["RNA"]]))
  
  cat("stored_filters_2\n")
  cat(str(stored_filters))
  cat("\n")
  
  
  ###### Remove high mito content cells (more or equal to 10% mt genes ------------------
  
  adata_sub_mito = subset(adata_sub_multiplet, subset = percent.mt < as.numeric(MITO_max))
  
  stored_filters['after_mito'] = length(colnames(adata_sub_mito[["RNA"]]))
  
  cat("stored_filters_3\n")
  cat(str(stored_filters))
  cat("\n")
  
  ### Redefine adata as the object with all the previous filters implemented ---------------
  
  adata <- adata_sub_mito
  
  ### Compute intermediate bulk metrics of filtered object ----------------
  
  DefaultAssay(adata) <- 'ATAC'
  adata <- TSSEnrichment(adata)
  stored_filters['median_TSSe'] = median(adata[[]][,'TSS.enrichment'])
  stored_filters['median_genesperCells_RNA'] = median(adata[[]][,'nFeature_RNA'])
  stored_filters['median_hq_atac_fragm'] = median(adata[[]][,'atac_fragments'])
  invisible(gc())
  
  
  #### Add sample name as the orig.ident field of the metadata ------
  
  adata@meta.data$orig.ident = sample_name
  
  
  metadata_adata<-adata[[]]
  
  cat("metadata_adata_REDEFINED\n")
  cat(str(metadata_adata))
  cat("\n")
  
 
  ########### PCA RNA analysis -------------------------
  
  cat("PCA RNA analysis\n")
  
  DefaultAssay(adata) <- 'RNA'
  suppressMessages(adata <- SCTransform(adata, verbose = FALSE) %>% RunPCA() %>% RunUMAP(dims=1:50, 
                                                                                         reduction.name='umap.rna', reduction.key='rnaUMAP_'))
  
  
  ########### PCA ATAC analysis -------------------------
  
  cat("PCA ATAC analysis\n")
  
  #We exclude the first dimension as this is typically correlated with sequencing depth
  DefaultAssay(adata) <- 'ATAC'
  adata <- RunTFIDF(adata)
  adata <- FindTopFeatures(adata, min.cutoff='q0')
  adata <- RunSVD(adata)
  adata <- RunUMAP(adata, reduction='lsi', dims=2:50, reduction.name='umap.atac', reduction.key='atacUMAP_')
  
  
  ################ Multimodal analysis using both RNA and ATAC layers for clusters ------------------------------------
  
  cat("Multimodal analysis using both RNA and ATAC layers for clusters\n")
  
  adata <- FindMultiModalNeighbors(adata, reduction.list=list('pca', 'lsi'), dims.list=list(1:50, 2:50))
  adata <- RunUMAP(adata, nn.name='weighted.nn', reduction.name='umap.wnn', reduction.key='wnnUMAP_')
  adata <- FindClusters(adata, graph.name='wsnn', algorithm=4, resolution = .5, verbose=FALSE)
  
   
  #### add Larry deconvolution results to metadata ------------------------------
  
  Uniquely_genotyped_larry_barcodes_assignments_sel<-Uniquely_genotyped_larry_barcodes_assignments[which(Uniquely_genotyped_larry_barcodes_assignments$sample == sample_name),]
  
  cat("Uniquely_genotyped_larry_barcodes_assignments_sel_0\n")
  cat(str(Uniquely_genotyped_larry_barcodes_assignments_sel))
  cat("\n")
  
  Uniquely_genotyped_larry_barcodes_assignments_sel_FILTERED<-Uniquely_genotyped_larry_barcodes_assignments_sel[which(Uniquely_genotyped_larry_barcodes_assignments_sel$Number_of_UMIS >= Threshold_UMIS_per_cell),]
  
  cat("Uniquely_genotyped_larry_barcodes_assignments_sel_FILTERED_0\n")
  cat(str(Uniquely_genotyped_larry_barcodes_assignments_sel_FILTERED))
  cat("\n")
  
  Uniquely_genotyped_larry_barcodes_assignments_sel_FILTERED$CellBC_adapted<-gsub("^[^:]+:[^:]+:","",Uniquely_genotyped_larry_barcodes_assignments_sel_FILTERED$CellBC)
  
  cat("Uniquely_genotyped_larry_barcodes_assignments_sel_FILTERED_1\n")
  cat(str(Uniquely_genotyped_larry_barcodes_assignments_sel_FILTERED))
  cat("\n")
  
  ### Add GFPbc_attribution
  
  gfp1      = setNames(Uniquely_genotyped_larry_barcodes_assignments_sel_FILTERED$GFPbc_attribution, 
                       Uniquely_genotyped_larry_barcodes_assignments_sel_FILTERED$CellBC_adapted)
  
  cat("gfp1_0\n")
  cat(str(gfp1))
  cat("\n")
  
  
  adata@meta.data$No_assigned_GFPbc = gfp1[rownames(adata@meta.data)]
  
  ### Add GFPbc
  
  gfp2      = setNames(as.character(Uniquely_genotyped_larry_barcodes_assignments_sel_FILTERED$GFPbc), 
                       Uniquely_genotyped_larry_barcodes_assignments_sel_FILTERED$CellBC_adapted)
  
  cat("gfp2_0\n")
  cat(str(gfp2))
  cat("\n")
  
  
  adata@meta.data$Assigned_GFPbc = gfp2[rownames(adata@meta.data)]
  

  metadata_adata<-adata[[]]
  
  cat("metadata_adata_REDEFINED_AFTER_GFPbc_assignation\n")
  cat(str(metadata_adata))
  cat("\n")
  

  adata@meta.data$Assigned_GFPgenotype<-NA
  
  adata@meta.data$Assigned_GFPgenotype[which(adata@meta.data$Assigned_GFPbc%in%c("chrGFP_WTA","chrGFP_WTB","chrGFP_WTC"))]<-'WT'
  adata@meta.data$Assigned_GFPgenotype[which(adata@meta.data$Assigned_GFPbc%in%c("chrGFP_rs1","chrGFP_rs2","chrGFP_rs3"))]<-'rsCHEK2'
  adata@meta.data$Assigned_GFPgenotype[which(adata@meta.data$Assigned_GFPbc%in%c("chrGFP_R882H1","chrGFP_R882H2","chrGFP_R882H3"))]<-'R882H'
  adata@meta.data$Assigned_GFPgenotype[which(adata@meta.data$Assigned_GFPbc%in%c("chrGFP_rs_R882H1","chrGFP_rs_R882H2","chrGFP_rs_R882H3"))]<-'Double_mutants'
  
  metadata_adata<-adata[[]]
  
  cat("metadata_adata_REDEFINED_AFTER_Assigned_GFPgenotype_assignation\n")
  cat(str(metadata_adata))
  cat("\n")
  
  adata@meta.data$Assigned_GFPbc<-factor(adata@meta.data$Assigned_GFPbc,
                                         levels = c("chrGFP_WTA","chrGFP_WTB","chrGFP_WTC",
                                                    "chrGFP_rs1","chrGFP_rs2","chrGFP_rs3",
                                                    "chrGFP_R882H1","chrGFP_R882H2","chrGFP_R882H3",
                                                    "chrGFP_rs_R882H1","chrGFP_rs_R882H2","chrGFP_rs_R882H3","No_GFPbcs"),
                                         ordered=T)
  
  
  
  adata@meta.data$Assigned_GFPgenotype<-factor(adata@meta.data$Assigned_GFPgenotype,
                                         levels = c("WT","rsCHEK2","R882H","Double_mutants"),
                                         ordered=T)
  
  
  metadata_adata<-adata[[]]
  
  cat("metadata_adata_REDEFINED_FINAL\n")
  cat(str(metadata_adata))
  cat("\n")
  
  #### adata_assigned -------
  
  adata_assigned = subset(adata, subset = No_assigned_GFPbc ==1 )
  
  
  #### SAVE object and selected barcodes --------------
  
  saveRDS(adata, file = file.path(output_dir,'preliminary_filtered.rds'))
  filtered_bcs <- colnames(adata[["RNA"]])
  write(filtered_bcs, file=(file.path(output_dir,'keep_barcodes_step1.txt')),sep='\n')
  
  
  #### output metrics ----------------------------------------------------------------------
  
  stored_filters["scDBL_RNA"] <- sum(adata@meta.data$scDblFinder.class == 'doublet')
  stored_filters["scDBL_ATAC"]<- sum(adata@meta.data$scDblFinder.class_atac == 'doublet')
  stored_filters["scDBL_both"]<- sum(adata@meta.data$DBL_comb == 'R:doublet:A:doublet')
  stored_filters["genotyped"] <- length(colnames(adata_assigned[["RNA"]]))
  
  write.table(data.frame(stored_filters), (file.path(output_dir,'barcodes_stats.tsv')),
              sep="\t", quote=F, col.names=FALSE)
  
}

graph_function = function(option_list)
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
  
  #### READ and transform master_path ----
  
  master_path = opt$master_path
  
  cat("master_path_\n")
  cat(sprintf(as.character(master_path)))
  cat("\n")
  
  sample_dir<-paste(master_path,sample_name,'/','outs','/',sep='')
  
  cat("sample_dir_\n")
  cat(sprintf(as.character(sample_dir)))
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
  
  path_processing_outputs = paste(out,sample_name,'/',sep='')
  
  if (file.exists(path_processing_outputs)){
    
    
  }else{
    
    dir.create(file.path(path_processing_outputs))
    
  }#path_processing_outputs
  
  
  output_dir = paste(path_processing_outputs,'intermediate','/',sep='')
  
  if (file.exists(output_dir)){
    
    
  }else{
    
    dir.create(file.path(output_dir))
    
  }#output_dir
  
  cat("output_dir\n")
  cat(sprintf(as.character(output_dir)))
  cat("\n")
  
  ##### Read the preprocessed object ------
  
  setwd(output_dir)
  
  adata<-readRDS(file="preliminary_filtered.rds")
  
  cat("adata_0\n")
  cat(str(adata))
  cat("\n")
  
  #### adata_assigned -------
  
  adata_assigned = subset(adata, subset = No_assigned_GFPbc ==1)
  
  ################# Barcode rank plots based on cellranger cell estimate ---------------------
  
  cat("Barcode rank plots based on cellranger cell estimate\n")
  
  sample_color = as.numeric(substring(sample_name, 7))
  ncells = read.csv(file.path(sample_dir, "summary.csv"))$Estimated.number.of.cells
  
  setwd(output_dir)
  
  png(file.path(output_dir,'Barcodes_rank_plots.png'), width = 1000, height = 500)
    par(mfrow=c(1,2), mar=c(6,6,6,6))
    bc_metrics <- read.csv(file.path(sample_dir, 'per_barcode_metrics.csv'),  , stringsAsFactors=1)
    plot( sort(bc_metrics$gex_umis_count,decreasing = T), log='xy', type='l', main=paste(sample_name,"RNA"),
          xlab="Barcode Rank", ylab="Number RNA UMIs", col="gray" ,lwd=3)
    lines( sort(bc_metrics$gex_umis_count,decreasing = T)[1:ncells],  col="blue4" ,lwd=3)
    grid()
    plot(sort(bc_metrics$atac_fragments, decreasing = T), log='xy', type='l', main=paste(sample_name,"ATAC"),
         xlab="Barcode Rank", ylab="Number ATAC fragment", col='gray', lwd=3)
    lines( sort(bc_metrics$atac_fragments,decreasing = T)[1:ncells],  col="dodgerblue" ,lwd=3)
    grid()
  dev.off()
  invisible(gc())
  
  ########### Plot Intermediate_metrics ---------------
  
  cat("Plot Intermediate_metrics\n")
  
  
  setwd(output_dir)
  
  Idents(adata) <- "orig.ident"
  
  # options(repr.plot.width = 10, repr.plot.height = 5)
  # png(file.path(output_dir,'Intermediate_metrics.png'), width =800, height = 400)
  #   VlnPlot(adata, features = c("nCount_ATAC", "nCount_RNA", "percent.mt",'TSS.enrichment'),  
  #           ncol = 4, cols=sample_color,
  #           log = TRUE, pt.size = 0)+ NoLegend()
  # dev.off()
  
  
 
  
  Vln_file<-VlnPlot(adata, features = c("nCount_ATAC", "nCount_RNA", "percent.mt",'TSS.enrichment'),  
          ncol = 4, cols=sample_color,
          log = TRUE, pt.size = 0)+ NoLegend()
  
  # saveRDS(Vln_file,file="Intermediate_metrics_graph.rds")
  
  
  graph <- Vln_file +
    plot_layout( guides = "collect") & theme(text = element_text(size = 8, family = "Arial"))

  # ggsave(p, filename = "plot.pdf", device=cairo_pdf, height=20, width=16, unit="cm", path=NULL, scale=1)
  # ggsave(p, filename = "plot.svg", device=svg, height=20, width=16, unit="cm", path=NULL, scale=1)

  ggsave(graph, filename = paste("Intermediate_metrics","_plot.pdf", sep=''), device=cairo_pdf, height=20, width=16, unit="cm", path=NULL, scale=1)
  
  
  # invisible(gc())
  
  ########### Plots of UMAPS and stats ---------------
  
  cat("Plots of UMAPS and stats\n")
  
  
  
  setwd(output_dir)
  
  options(repr.plot.width = 12, repr.plot.height = 5)
  png(file.path(output_dir,'Intermediate_UMAPs_clusters.png'), width =1000, height = 400)
    p1 <- DimPlot(adata, reduction = "umap.rna", group.by = "seurat_clusters", label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle("RNA")
    p2 <- DimPlot(adata, reduction = "umap.atac", group.by = "seurat_clusters", label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle("ATAC")
    p3 <- DimPlot(adata, reduction = "umap.wnn", group.by = "seurat_clusters", label = TRUE, label.size = 2.5, repel = TRUE) + ggtitle("WNN")
    p1 + p2 + p3 & NoLegend() & theme(plot.title = element_text(hjust = 0.5))
  dev.off()
  
  graph<-p1 + p2 + p3 & NoLegend() & theme(plot.title = element_text(hjust = 0.5))
  
  saveRDS(graph,file="UMAPS_and_stats_graph.rds")
  
 
  
  
  ########### Plot Intermediate_metrics and marker genes ---------------
  
  cat("Plot Intermediate_metrics and marker genes\n")
  
  

  setwd(output_dir)
  
  options(repr.plot.width = 14, repr.plot.height = 8)
  png(file.path(output_dir,'Intermediate_UMAPs_qc.png'), width =1500, height = 1800)
    p5 <- FeaturePlot(adata, features = c("NANOG",'SOX4', 'EOMES', 'TWIST1'),
                      reduction = 'umap.wnn', 
                      cols = c("lightgrey","darkgreen"), ncol = 4)
    p5_1 <- FeaturePlot(adata, features = c('GYPA','HBA2','HBZ','ITGA2B'),
                        reduction = 'umap.wnn', 
                        cols = c("lightgrey","darkgreen"), ncol = 4)
    p5_2 <- FeaturePlot(adata, features = c('GP1BA','KIF15','STIL','ANAPC7'),
                        reduction = 'umap.wnn', 
                        cols = c("lightgrey","darkgreen"), ncol = 4)
    p6 <- FeaturePlot(adata, features = c("nCount_SCT", "nCount_RNA", "nCount_ATAC",'TSS.enrichment'), ncol = 4,
                      reduction = 'umap.wnn')
    p7 <- FeaturePlot(adata, features = c("nFeature_SCT", "nFeature_RNA", "nFeature_ATAC",'percent.mt'), ncol = 4,
                      reduction = 'umap.wnn')
    p6 / p7 / p5 / p5_1 /p5_2
  dev.off()
  
  graph<-p6 / p7 / p5 / p5_1 /p5_2
  
  saveRDS(graph,file="marker_genes_graph.rds")
  
 
  ### Plot Intermediate_UMAPs_doublets.png
  
  cat("Plot Intermediate_UMAPs_doublets.png\n")
  
  
  setwd(output_dir)
  
  png(file.path(output_dir,'Intermediate_UMAPs_doublets.png'), width =1000, height = 300)
  
    p1<- DimPlot(adata, reduction = "umap.rna", group.by = "scDblFinder.class", label.size = 2.5, repel = TRUE) + ggtitle("RNA")+ scale_color_manual(values = c("gray", "black"))
    p2<- DimPlot(adata, reduction = "umap.atac", group.by = "scDblFinder.class_atac", label.size = 2.5, repel = TRUE) + ggtitle("ATAC") + scale_color_manual(values = c("gray", "black"))
    p3<-DimPlot(adata, reduction = "umap.wnn", group.by = "DBL_comb", label.size = 2.5, repel = TRUE) + ggtitle("WNN") + scale_color_manual(values = c("black", "orange", "gold","gray"))
    
    p1 + p2 +p3
  dev.off()
  
  graph<-p1 + p2 +p3
  
  
  saveRDS(graph,file="Intermediate_UMAPs_doublets_graph.rds")
  
  #### Plot Intermediate_UMAPs_unique_GFP_bc
  
  cat("Plot Intermediate_UMAPs_unique_GFP_bc\n")
  
  
  setwd(output_dir)
  
  png(file.path(output_dir,'Intermediate_UMAPs_unique_GFP_bc.png'), width =1000, height = 600)
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
  
  #### READ and transform master_path ----
  
  master_path = opt$master_path
  
  cat("master_path_\n")
  cat(sprintf(as.character(master_path)))
  cat("\n")
  
  sample_dir<-paste(master_path,sample_name,'/','outs','/',sep='')
  
  cat("sample_dir_\n")
  cat(sprintf(as.character(sample_dir)))
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
  
  path_processing_outputs = paste(out,sample_name,'/',sep='')
  
  if (file.exists(path_processing_outputs)){
    
    
  }else{
    
    dir.create(file.path(path_processing_outputs))
    
  }#path_processing_outputs
  
  
  output_dir = paste(path_processing_outputs,'intermediate','/',sep='')
  
  if (file.exists(output_dir)){
    
    
  }else{
    
    dir.create(file.path(output_dir))
    
  }#output_dir
  
  cat("output_dir\n")
  cat(sprintf(as.character(output_dir)))
  cat("\n")
  
  ##### Read the preprocessed object ------
  
  setwd(output_dir)
  
  
  
  ############## Retrieve graphs  -----------------
  
  
  file_list <- list.files(path=output_dir, include.dirs = FALSE)
  
  
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
    
    setwd(output_dir)
    
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
    
    
    setwd(output_dir)
    
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
  
  # setwd(output_dir)
  # 
  # p <- ((List_RESULTS[[1]] / List_RESULTS[[2]] / List_RESULTS[[3]]/ List_RESULTS[[4]]/ List_RESULTS[[5]])) +
  #   plot_annotation(tag_levels = "a") +
  #   plot_layout( guides = "collect") & theme(text = element_text(size = 8, family = "Arial"))
  # 
  # ggsave(p, filename = "plot.pdf", device=cairo_pdf, height=20, width=16, unit="cm", path=NULL, scale=1)
  # ggsave(p, filename = "plot.svg", device=svg, height=20, width=16, unit="cm", path=NULL, scale=1)
  
  
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
    make_option(c("--Uniquely_genotyped_larry_barcodes_assignments"), type="character", default=NULL, 
                metavar="type", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--Threshold_UMIS_per_cell"), type="numeric", default=NULL, 
                metavar="type", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--sample_name"), type="character", default=NULL, 
                metavar="type", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--master_path"), type="character", default=NULL, 
                metavar="type", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--rna_min_features"), type="numeric", default=NULL, 
                metavar="type", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--atac_min_fragments"), type="numeric", default=NULL, 
                metavar="type", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--MITO_max"), type="numeric", default=NULL, 
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
  
  create_the_preliminary_Seurat_object(opt)
  graph_function(opt)
  patchwork_function(opt)
 

}


###########################################################################

system.time( main() )
  