


.libPaths()

assign(".lib.loc", "/home/manuel.tardaguila/conda_envs/multiome_QC_DEF/lib/R/library", envir = environment(.libPaths))

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
suppressMessages(library(SeuratData))
suppressMessages(library(SeuratDisk))



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
  
 
  #### Read filtered object by doublets -----
  
  
  adata<-readRDS(file=opt$Seurat_object)
  
  cat("adata_0\n")
  # cat(str(adata))
  cat("\n")
  
   
   
   
   adata@meta.data$cell_type<-as.character(adata@meta.data$refined_annotation_majority_vote)
   adata@meta.data$clone_line<-as.character(adata@meta.data$clone_line)
   adata@meta.data$Genotype<-as.character(adata@meta.data$Genotype)

   met<-adata[[]]
   
   cat("met_0\n")
   cat(str(met))
   cat("\n")

  ##### Reduce the Seurat object to h5ad with RNA counts not corrected by CellBender It doesn't work with CellBender corrected counts------

  DefaultAssay(adata)<-'RNA_raw'
  RNA_only<-DietSeurat(adata, assays = "RNA_raw")


  setwd(out)
  
  unlink(c("RNA.h5Seurat","RNA.h5ad"))

  SaveH5Seurat(RNA_only, filename = "RNA.h5Seurat")
  Convert("RNA.h5Seurat", dest = "h5ad")


  ##### Reduce the Seurat object to h5ad with ATAC counts not corrected by CellBender It doesn't work with CellBender corrected counts------

  DefaultAssay(adata)<-'ATAC_by_refined_annotation'
  
  
  ATAC_only<-DietSeurat(adata, assays = "ATAC_by_refined_annotation")
  
  
  Peaks<-Features(adata)
  
  cat("Peaks_0\n")
  cat(str(Peaks))
  cat("\n")
  
  Peaks<-unique(Peaks)
  
  cat("Peaks_1\n")
  cat(str(Peaks))
  cat("\n")
  
  
  tmp.gather<- data.frame(matrix(vector(), length(Peaks), 4,
                                 dimnames=list(c(),
                                               c("chr","start","end","name"))),
                          stringsAsFactors=F)
  

  tmp.gather$name<-Peaks
  tmp.gather$chr<-gsub("-.+$","",Peaks)
  
  tmp.gather$start<-gsub("^[^-]+-","",Peaks)
  tmp.gather$start<-as.integer(gsub("-.+$","",tmp.gather$start))
  tmp.gather$end<-as.integer(gsub("^[^-]+-[^-]+-","",Peaks))
  
  
  cat("tmp.gather_0\n")
  cat(str(tmp.gather))
  cat("\n")
  
  tmp.gather$chr<-factor(tmp.gather$chr,
                           levels=c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11",
                                    "chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20","chr21",
                                    "chr22","chr23","chrX","chrY"), ordered=T)
  
  tmp.gather<-tmp.gather[order(tmp.gather$chr, tmp.gather$start, decreasing = F),]
  
  cat("tmp.gather_1\n")
  cat(str(tmp.gather))
  cat("\n")
  

  # Peaks$chr<-gsub("-.+$","",)


  setwd(out)
  
  unlink(c("ATAC.h5Seurat","ATAC.h5ad"))

  SaveH5Seurat(ATAC_only, filename = "ATAC.h5Seurat")
  Convert("ATAC.h5Seurat", dest = "h5ad")

  
  write.table(tmp.gather, file="Peaks.bed", sep="\t", row.names = F,quote = F,col.names = F)
  
  
  
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
    make_option(c("--Seurat_object"), type="character", default=NULL, 
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
  