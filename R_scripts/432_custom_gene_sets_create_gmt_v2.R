
suppressMessages(library("plyr", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("data.table", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("crayon", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("withr", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("ggplot2", lib.loc = "/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("farver", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("labeling", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("optparse", lib.loc = "/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("dplyr", lib.loc = "/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("withr", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("backports", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("broom", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("rstudioapi", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("cli", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("tzdb", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("svglite", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("ggeasy", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("sandwich", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("digest", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("tidyverse", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("RColorBrewer", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("svglite", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("cowplot", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("org.Hs.eg.db", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("ActivePathways", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("splitstackshape", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("dorothea", lib="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("decoupleR", lib="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("rWikiPathways", lib="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))

suppressMessages(library("splitstackshape", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))


opt = NULL

options(warn = 1)

create_gmt = function(option_list)
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
  
  cat("OUT_\n")
  cat(sprintf(as.character(out)))
  cat("\n")
  
  
  #### Table_of_gene_sets_Manual_curation ----
  
  Table_of_gene_sets<-read.table(opt$Table_of_gene_sets, sep="\t", header=T)
  

  cat("Table_of_gene_sets_0\n")
  cat(str(Table_of_gene_sets))
  cat("\n")
  cat(str(unique(Table_of_gene_sets$GeneSet)))
  cat("\n")
  
  Table_of_gene_sets_long<-unique(as.data.frame(cSplit(Table_of_gene_sets,sep = ',', direction = "long",
                                               splitCols = "Genes"),stringsAsFactors=F))
  
  cat("Table_of_gene_sets_long_0\n")
  cat(str(Table_of_gene_sets_long))
  cat("\n")
  cat(str(unique(Table_of_gene_sets_long$GeneSet)))
  cat("\n")
  
  
  
 
  
  ### Retrieve the ENTREZID and ensembl_gene_id ----
  
  
  Table_of_gene_sets_long$ENTREZID<-mapIds(org.Hs.eg.db, keys=Table_of_gene_sets_long$Genes, keytype="SYMBOL",column="ENTREZID")
  Table_of_gene_sets_long$ensembl_gene_id<-mapIds(org.Hs.eg.db, keys=Table_of_gene_sets_long$Genes, keytype="SYMBOL",column="ENSEMBL")
  
  cat("Table_of_gene_sets_long_1\n")
  cat(str(Table_of_gene_sets_long))
  cat("\n")
 
  Table_of_gene_sets_long_NO_NA<-Table_of_gene_sets_long[!is.na(Table_of_gene_sets_long$ENTREZID),]
  
  cat("Table_of_gene_sets_long_NO_NA_0\n")
  cat(str(Table_of_gene_sets_long_NO_NA))
  cat("\n")
  
  #### Prepare the file for gmt export ----
  
  for_gmt_Table_of_gene_sets_long_NO_NA<-unique(Table_of_gene_sets_long_NO_NA[,c(which(colnames(Table_of_gene_sets_long_NO_NA) == 'GeneSet'),which(colnames(Table_of_gene_sets_long_NO_NA) == 'ENTREZID'))])
  
  colnames(for_gmt_Table_of_gene_sets_long_NO_NA)[which(colnames(for_gmt_Table_of_gene_sets_long_NO_NA) == 'GeneSet')]<-'id'
  colnames(for_gmt_Table_of_gene_sets_long_NO_NA)[which(colnames(for_gmt_Table_of_gene_sets_long_NO_NA) == 'ENTREZID')]<-'gene'
  
  cat("for_gmt_Table_of_gene_sets_long_NO_NA_0\n")
  cat(str(for_gmt_Table_of_gene_sets_long_NO_NA))
  cat("\n")
  
  
  for_gmt_Table_of_gene_sets_long_NO_NA$name<-for_gmt_Table_of_gene_sets_long_NO_NA$id
  
  
  
  cat("for_gmt_Table_of_gene_sets_long_NO_NA_1\n")
  cat(str(for_gmt_Table_of_gene_sets_long_NO_NA))
  cat("\n")
  cat(sprintf(as.character(names(summary(as.factor(for_gmt_Table_of_gene_sets_long_NO_NA$name))))))
  cat("\n")
  cat(sprintf(as.character(summary(as.factor(for_gmt_Table_of_gene_sets_long_NO_NA$name)))))
  cat("\n")
  
  
  #### Reorder
  
  indx.reorder<-c(which(colnames(for_gmt_Table_of_gene_sets_long_NO_NA) == 'id'),which(colnames(for_gmt_Table_of_gene_sets_long_NO_NA) == 'name'),which(colnames(for_gmt_Table_of_gene_sets_long_NO_NA) == 'gene'))
  
  
  for_gmt_Table_of_gene_sets_long_NO_NA_reordered<-for_gmt_Table_of_gene_sets_long_NO_NA[,indx.reorder]
  
  cat("for_gmt_Table_of_gene_sets_long_NO_NA_reordered_0\n")
  cat(str(for_gmt_Table_of_gene_sets_long_NO_NA_reordered))
  cat("\n")
 
 ### save as gmt ----
  
  setwd(out)
  
  writeGMT(for_gmt_Table_of_gene_sets_long_NO_NA_reordered, paste("Custom_Soranzo",'_Hs.entrez.gmt',sep=''))
  
  
  
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
    make_option(c("--Table_of_gene_sets"), type="character", default=NULL, 
                metavar="type", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--type"), type="character", default=NULL, 
                metavar="type", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--out"), type="character", default=NULL, 
                metavar="filename", 
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
  
  create_gmt(opt)

  
}


###########################################################################

system.time( main() )