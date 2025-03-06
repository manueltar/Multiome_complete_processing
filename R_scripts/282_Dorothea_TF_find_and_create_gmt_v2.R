
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



opt = NULL

options(warn = 1)

Dorothea_and_gmt = function(option_list)
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
  
  #### READ and transform Dorothea_confidence ----
  
  Dorothea_confidence = unlist(strsplit(opt$Dorothea_confidence, split="_"))
  
  cat("Dorothea_confidence_\n")
  cat(sprintf(as.character(Dorothea_confidence)))
  cat("\n")
  
  #### READ and transform TF_annotation ----
  
  TF_annotation = unlist(strsplit(opt$TF_annotation, split=","))
  
  cat("TF_annotation_\n")
  cat(sprintf(as.character(TF_annotation)))
  cat("\n")
  
  #### READ and transform selected_TFs ----
  
  selected_TFs = unlist(strsplit(opt$selected_TFs, split=","))
  
  cat("selected_TFs_\n")
  cat(sprintf(as.character(selected_TFs)))
  cat("\n")
  
  #### Table_S6_Manual_curation ----
  
  Table_S6<-readRDS(opt$Table_S6)
  
  Table_S6_subset<-droplevels(Table_S6[which(Table_S6$Manual_curation == 'R in candidate'),])
  
  cat("Table_S6_subset_0\n")
  cat(str(Table_S6_subset))
  cat("\n")
  cat(str(unique(Table_S6_subset$VAR)))
  cat("\n")
  
  Table_S6_subset$TF_CLASS<-NA
  
  
  Table_S6_subset$TF_CLASS[which(Table_S6_subset$VAR%in%TF_annotation)]<-'TF'
  Table_S6_subset$TF_CLASS[-which(Table_S6_subset$VAR%in%TF_annotation)]<-'NO_TF'
  
  Table_S6_subset$TF_CLASS<-factor(Table_S6_subset$TF_CLASS, levels=c('NO_TF','TF'))
  
  
  cat("Table_S6_subset_1\n")
  cat(str(Table_S6_subset))
  cat("\n")
  cat(str(unique(Table_S6_subset$VAR)))
  cat("\n")
  cat(sprintf(as.character(names(summary(Table_S6_subset$TF_CLASS)))))
  cat("\n")
  cat(sprintf(as.character(summary(Table_S6_subset$TF_CLASS))))
  cat("\n")
  
  #### TF genes ----
  
  Table_S6_TF<-droplevels(Table_S6_subset[which(Table_S6_subset$TF_CLASS == 'TF'),])
  
  cat("Table_S6_TF_0\n")
  cat(str(Table_S6_TF))
  cat("\n")
  cat(str(unique(Table_S6_TF$VAR)))
  cat("\n")
  cat(sprintf(as.character(names(summary(Table_S6_TF$TF_CLASS)))))
  cat("\n")
  cat(sprintf(as.character(summary(Table_S6_TF$TF_CLASS))))
  cat("\n")
  
  
  TF_genes<-unique(c(unlist(strsplit(Table_S6_TF$Whole_blood_DE_HGNC_string, split=";")),
                     unlist(strsplit(Table_S6_TF$Monocyte_DE_HGNC_string, split=";")),
                     unlist(strsplit(Table_S6_TF$Tcell_DE_HGNC_string, split=";")),
                     unlist(strsplit(Table_S6_TF$Neutrophil_DE_HGNC_string, split=";")),
                     unlist(strsplit(Table_S6_TF$Whole_blood_DTU_HGNC_string, split=";")),
                     unlist(strsplit(Table_S6_TF$Monocyte_DTU_HGNC_string, split=";")),
                     unlist(strsplit(Table_S6_TF$Tcell_DTU_HGNC_string, split=";")),
                     unlist(strsplit(Table_S6_TF$Neutrophil_DTU_HGNC_string, split=";"))))
  
  cat("TF_genes_0\n")
  cat(str(TF_genes))
  cat("\n")
  cat(sprintf(as.character(TF_genes)))
  cat("\n")
  
  TF_genes_NO_NA<-TF_genes[!is.na(TF_genes)]
  
  
  cat("TF_genes_NO_NA_0\n")
  cat(str(TF_genes_NO_NA))
  cat("\n")
  cat(sprintf(as.character(TF_genes_NO_NA)))
  cat("\n")
  
  #### Extract Dorothea genes ----
  
  net <- decoupleR::get_dorothea(levels = c(Dorothea_confidence))
  
  
  cat("net_0\n")
  cat(str(net))
  cat("\n")
  
  net.df<-as.data.frame(net)
  
  cat("net.df_0\n")
  cat(str(net.df))
  cat("\n")
  
  
  #### indx grep the TF genes ----
  
  indx.int<-grep(paste(unique(c(TF_genes_NO_NA,selected_TFs)), collapse='|'), net.df$source)
  
  cat("indx.int_0\n")
  cat(str(indx.int))
  cat("\n")
  
  check_sources<-unique(net.df$source[indx.int])
  
  cat("check_sources\n")
  cat(sprintf(as.character(check_sources)))
  cat("\n")
  
  
  
  net.df_oI<-net.df[which(net.df$source%in%unique(c(TF_genes_NO_NA,selected_TFs))),]
  
  cat("net.df_oI_0\n")
  cat(str(net.df_oI))
  cat("\n")
  cat(sprintf(as.character(names(summary(as.factor(net.df_oI$source))))))
  cat("\n")
  cat(sprintf(as.character(summary(as.factor(net.df_oI$source)))))
  cat("\n")
  
  ### Retrieve the ENTREZID and ensembl_gene_id ----
  
  
  net.df_oI$ENTREZID<-mapIds(org.Hs.eg.db, keys=net.df_oI$target, keytype="SYMBOL",column="ENTREZID")
  net.df_oI$ensembl_gene_id<-mapIds(org.Hs.eg.db, keys=net.df_oI$target, keytype="SYMBOL",column="ENSEMBL")
  
  cat("net.df_oI_1\n")
  cat(str(net.df_oI))
  cat("\n")
 
  net.df_oI_NO_NA<-net.df_oI[!is.na(net.df_oI$ENTREZID),]
  
  cat("net.df_oI_NO_NA_0\n")
  cat(str(net.df_oI_NO_NA))
  cat("\n")
  
  #### Prepare the file for gmt export ----
  
  for_gmt_net.df_oI_NO_NA<-unique(net.df_oI_NO_NA[,c(which(colnames(net.df_oI_NO_NA) == 'source'),which(colnames(net.df_oI_NO_NA) == 'ENTREZID'))])
  
  colnames(for_gmt_net.df_oI_NO_NA)[which(colnames(for_gmt_net.df_oI_NO_NA) == 'source')]<-'id'
  colnames(for_gmt_net.df_oI_NO_NA)[which(colnames(for_gmt_net.df_oI_NO_NA) == 'ENTREZID')]<-'gene'
  
  cat("for_gmt_net.df_oI_NO_NA_0\n")
  cat(str(for_gmt_net.df_oI_NO_NA))
  cat("\n")
  
  
  
  Symbol_array<-unique(for_gmt_net.df_oI_NO_NA$id)
  
  
  for_gmt_net.df_oI_NO_NA$name<-NA
  
  DEBUG<-1
  
  for(i in 1:length(Symbol_array)){
    
    Symbol_array_sel<-Symbol_array[i]
    
    cat("--------------------------------------->\t")
    cat(sprintf(as.character(Symbol_array_sel)))
    cat("\n")
    
    indx.Symbol<-which(for_gmt_net.df_oI_NO_NA$id == Symbol_array_sel)
    
    if(DEBUG ==1)
    {
      cat("indx.Symbol_0\n")
      cat(str(indx.Symbol))
      cat("\n")
    }
    
    for_gmt_net.df_oI_NO_NA$name[indx.Symbol]<-paste('Dorothea_',paste(Dorothea_confidence, collapse=""),Symbol_array_sel,'targets', sep="_")
    for_gmt_net.df_oI_NO_NA$id[indx.Symbol]<-paste('Dorothea_',paste(Dorothea_confidence, collapse=""),Symbol_array_sel,'targets', sep="_")
    
    
    
    
  }#i in 1:length(Symbol_array)
  
  
  cat("for_gmt_net.df_oI_NO_NA_1\n")
  cat(str(for_gmt_net.df_oI_NO_NA))
  cat("\n")
  cat(sprintf(as.character(names(summary(as.factor(for_gmt_net.df_oI_NO_NA$name))))))
  cat("\n")
  cat(sprintf(as.character(summary(as.factor(for_gmt_net.df_oI_NO_NA$name)))))
  cat("\n")
  
  
  #### Reorder
  
  indx.reorder<-c(which(colnames(for_gmt_net.df_oI_NO_NA) == 'id'),which(colnames(for_gmt_net.df_oI_NO_NA) == 'name'),which(colnames(for_gmt_net.df_oI_NO_NA) == 'gene'))
  
  
  for_gmt_net.df_oI_NO_NA_reordered<-for_gmt_net.df_oI_NO_NA[,indx.reorder]
  
  cat("for_gmt_net.df_oI_NO_NA_reordered_0\n")
  cat(str(for_gmt_net.df_oI_NO_NA_reordered))
  cat("\n")
 
 ### save as gmt ----
  
  setwd(out)
  
  writeGMT(for_gmt_net.df_oI_NO_NA_reordered, paste("Dorothea_",paste(Dorothea_confidence, collapse=""),'_',paste(Symbol_array, collapse="_"),'_Hs.entrez.gmt',sep=''))
  
  
  
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
    make_option(c("--Table_S6"), type="character", default=NULL, 
                metavar="type", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--selected_TFs"), type="character", default=NULL, 
                metavar="filename", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--Dorothea_confidence"), type="character", default=NULL, 
                metavar="filename", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--TF_annotation"), type="character", default=NULL, 
                metavar="filename", 
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
  
  Dorothea_and_gmt(opt)

  
}


###########################################################################

system.time( main() )