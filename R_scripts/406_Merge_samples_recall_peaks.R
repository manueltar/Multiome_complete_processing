
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




opt = NULL

options(warn = -1)

merge_and_recall_peaks = function(option_list)
{
  
  opt_in = option_list
  opt <<- option_list
  
  cat("All options:\n")
  printList(opt)
  
  
  #### READ and transform sample_array ----
  
  sample_array = unlist(strsplit(opt$sample_array, split=","))
  
  cat("sample_array_\n")
  cat(sprintf(as.character(sample_array)))
  cat("\n")
  
  #### READ and transform type ----
  
  type = opt$type
  
  cat("TYPE_\n")
  cat(sprintf(as.character(type)))
  cat("\n")
  
  
  #### LOOP TO READ in the pre merge object per sample and transform out ----
  
  out = opt$out
  
  cat("out_\n")
  cat(sprintf(as.character(out)))
  cat("\n")
  
  adatas <- list()
  
  for(i in 1:length(sample_array)){
    
    sample_array_sel<-sample_array[i]
    
    cat("---------------------------------------->\t") 
    cat(sprintf(as.character(sample_array_sel)))
    cat("\n")    
    
    path_processing_outputs = paste(out,sample_array_sel,'/',sep='')
    
    if (file.exists(path_processing_outputs)){
      
      
    }else{
      
      dir.create(file.path(path_processing_outputs))
      
    }#path_processing_outputs
    
    
    pre_merge_dir = paste(path_processing_outputs,'pre_merge','/',sep='')
    
    if (file.exists(pre_merge_dir)){
      
      
    }else{
      
      dir.create(file.path(pre_merge_dir))
      
    }#pre_merge_dir
    
    setwd(pre_merge_dir)
    
    adata <- readRDS("pre_merged.rds")
    
    adata@meta.data$orig.ident = sample_array_sel
    adatas[[sample_array_sel]] <- adata
    DefaultAssay(adatas[[sample_array_sel]]) <- "RNA"
    adatas[[sample_array_sel]][['ATAC']] <- NULL
 
    
    
    snATAC_matrices_dir = paste(path_processing_outputs,'snATAC_matrices','/',sep='')
    
    if (file.exists(snATAC_matrices_dir)){
      
      
    }else{
      
      dir.create(file.path(snATAC_matrices_dir))
      
    }#snATAC_matrices_dir
    
    setwd(snATAC_matrices_dir)
    
    
    filename<-paste(sample_array_sel, '_snATAC_pipeline_job.long_fmt_mtx.txt.gz', sep='')
    
    cat("filename_\n")
    cat(sprintf(as.character(filename)))
    cat("\n")
    
    lfmat      = read.table(filename)
    lfmat$V1   = paste(sample_array_sel,lfmat$V1, sep="_" )    
    if(sample_array_sel==sample_array[1]){
      LFM = lfmat
      
      }else{
        
        # Do nothing
        
      }# sample_array_sel==sample_array[1]
    
    LFM = rbind(LFM, lfmat)   
    
  }#i in 1:length(sample_array)
  
  
  merged = merge(x =adatas[[1]], y=adatas[2:4], add.cell.ids = sample_array )
  
  cat("merged_0\n")
  cat(str(merged))
  cat("\n")
  
  merged[["RNA"]] <-JoinLayers(merged[["RNA"]])
  
  cat("merged_1\n")
  cat(str(merged))
  cat("\n")
  
  rm(adatas)
  
  gc()
  
  atac_sm <- with(LFM,
                  sparseMatrix(i=as.numeric(as.factor(V2)), j=as.numeric(as.factor(V1)), 
                               x=V3, dimnames=list(levels(as.factor(V2)), levels(as.factor(V1)))))
  
  cat("atac_sm_1\n")
  cat(str(atac_sm))
  cat("\n")
  
  
 
  ############################################################
  #create the new chromatin assay object and add to Seurat object
  ############################################################
  
  atac_sm       <- atac_sm[,colnames(merged)]
  grange.counts <- StringToGRanges(rownames(atac_sm), sep = c(':', '-'))
  grange.use    <- seqnames(grange.counts) %in% standardChromosomes(grange.counts)
  atac_sm       <- atac_sm[as.vector(grange.use), ]
  suppressMessages(annotations <- GetGRangesFromEnsDb(ensdb=EnsDb.Hsapiens.v86))
  seqlevelsStyle(annotations)  <- 'UCSC'
  genome(annotations)          <- 'hg38'
  
  suppressWarnings(chrom_assay <- CreateChromatinAssay(counts=atac_sm, sep=c(':', '-'), 
                                                       genome='hg38', fragments=opt$frag_file, 
                                                       min.cells=-1, min.features=-1, 
                                                       annotation=annotations))
  
  
  
  merged[['ATAC']] <- chrom_assay
  
  cat("merged_2\n")
  cat(str(merged))
  cat("\n")
  
  ###### SAVE -----
  
  setwd(out)
  
  saveRDS(merged, file = 'merged_unprocessed.rds')
  
  
  
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
    make_option(c("--sample_array"), type="character", default=NULL, 
                metavar="type", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--frag_file"), type="character", default=NULL, 
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
  