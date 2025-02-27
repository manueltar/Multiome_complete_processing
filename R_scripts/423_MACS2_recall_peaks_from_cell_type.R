
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

MACS2_call_peaks = function(option_list)
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
  
  #### Read seurat_object_annotated -----
  
  
  adata<-readRDS(file=opt$seurat_object_annotated)
  
  # cat("adata_0\n")
  # cat(str(adata))
  # cat("\n")
  
  
  
  ## Call Peaks and make new peak matrix
  
  
  DefaultAssay(adata) <- 'ATAC'
  peaks <- CallPeaks(
    object = adata,
    group.by = "refined_annotation_majority_vote",    
    macs2.path = "/group/soranzo/conda_envs/Manuel_macs2/bin/macs2")
  
  frag_file<-opt$frag_file
  
  Fragmobj <- CreateFragmentObject(frag_file,cells =Cells(adata))
  
  
  peakmat = FeatureMatrix(fragments = Fragmobj, features = peaks, cells = Cells(adata), process_n = 10000,
                          sep = c(":", "-"), verbose = TRUE)
  
  norm_chr = rownames(peakmat)[stringr::str_split_fixed(rownames(peakmat), ":",2)[,1] %in% 
                                 paste0("chr", c(1:22, "X", "Y"))]
  
  peakmat=peakmat[norm_chr,]
  
  suppressMessages(annotations <- GetGRangesFromEnsDb(ensdb=EnsDb.Hsapiens.v86))
  seqlevelsStyle(annotations)  <- 'UCSC'
  genome(annotations)          <- 'hg38'
  
  suppressWarnings(chrom_assay <- CreateChromatinAssay(counts=peakmat, sep=c(':', '-'), 
                                                       genome='hg38', fragments=Fragmobj, 
                                                       min.cells=-1, min.features=-1, 
                                                       annotation=annotations))
  
  
  
  adata[['ATAC_by_refined_annotation']] <- chrom_assay
  
  #### Save third clustered filtered object  ------------------
  
  
  setwd(out)
  
  saveRDS(adata,file="merged_clusters_after_genotyping_after_refined_annotation_new_peaks.rds")
  
  
  
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
    make_option(c("--seurat_object_annotated"), type="character", default=NULL, 
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
  
  MACS2_call_peaks(opt)

 

}


###########################################################################

system.time( main() )
  