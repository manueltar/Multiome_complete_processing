.libPaths()
.libPaths(new = c("/home/manuel.tardaguila/conda_envs/multiome_NEW_downstream_analysis/lib/R/library"))
.libPaths()
# sessionInfo()


suppressMessages(library(dplyr)) 
suppressMessages(library(ggplot2)) 
suppressMessages(library(Matrix)) 
suppressMessages(library(data.table)) 
suppressMessages(library(ggpubr)) 
suppressMessages(library(ggplot2))
suppressMessages(library(pheatmap))
suppressMessages(library("cowplot"))
suppressMessages(library("RColorBrewer"))
suppressMessages(library("plyr"))
suppressMessages(library("forcats"))
suppressMessages(library('ggeasy'))
suppressMessages(library('dplyr'))
suppressMessages(library("svglite"))
suppressMessages(library("ape"))
suppressMessages(library("ggforce"))
suppressMessages(library("tidyr"))
suppressMessages(library("tibble")) 
library("ggrepel")

library("optparse")
suppressMessages(library("splitstackshape")) 
suppressMessages(library("ggupset"))



opt = NULL

options(warn = 1)



multiVals <- function(x) paste(x,collapse=";")


line_plot_function = function(option_list)
{
  
  opt_in = option_list
  opt <<- option_list
  
  cat("All options:\n")
  printList(opt)
  
  
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
  
  #### READ and transform reference ----
  
  reference = opt$reference
  
  cat("reference_\n")
  cat(sprintf(as.character(reference)))
  cat("\n")
 
  #### READ selected_clone_lines ----
  
  selected_clone_lines = unlist(strsplit(opt$selected_clone_lines, split=","))
  
  cat("selected_clone_lines_0\n")
  cat(sprintf(as.character(selected_clone_lines)))
  cat("\n")
  

  #### READ and transform adata ----
  
  
  adata<-readRDS(file=opt$SeuratObject)
  
  cat("adata_0\n")
  # cat(str(adata))
  cat("\n")
  
  
  met<-adata[[]]
  
  cat("met_0\n")
  cat(str(met))
  cat("\n")
  
  
  met_sel<-met[which(met$clone_line%in%selected_clone_lines),]
  
  
  cat("met_sel_0\n")
  cat(str(met_sel))
  cat("\n")
  
  
  met_sel.dt<-data.table(met_sel,key=c("time_point","clone_line","Genotype"))
  

  Freq.table<-as.data.frame(met_sel.dt[,.(Freq=.N),by=key(met_sel.dt)], stringsAsFactors=F)
  
  cat("Freq.table_0\n")
  cat(str(Freq.table))
  cat("\n")
  
  
  indx.ref<-which(Freq.table$Genotype%in%reference)
  
  cat("indx.ref_1\n")
  cat(str(indx.ref))
  cat("\n")
  
  
  
  reference_df<-Freq.table[indx.ref,]
  
  cat("reference_df\n")
  str(reference_df)
  cat("\n")
  
  reference_df.dt<-data.table(reference_df,
                            key=c('time_point'))
  
  
  
  Mean_counts_relative_to_ref<-as.data.frame(reference_df.dt[,.(mean_cells_ref=mean(Freq)), by=key(reference_df.dt)], stringsAsFactors=F)
  
  
  cat("Mean_counts_relative_to_ref_0\n")
  str(Mean_counts_relative_to_ref)
  cat("\n")
  

  Freq.table<-merge(Freq.table,
                    Mean_counts_relative_to_ref,
                    by=c('time_point'))
  
  cat("Freq.table_AFTER_merge_with_reference\n")
  cat(str(Freq.table))
  cat("\n")
  
  Freq.table.dt<-data.table(Freq.table,
                            key=c('Genotype','time_point'))
  
  
  Mean_counts_relative<-as.data.frame(Freq.table.dt[,.(cells_relative_to_ref=mean(Freq/mean_cells_ref),
                                                       sd=sd(Freq/mean_cells_ref),
                                                       samples=.N), by=key(Freq.table.dt)], stringsAsFactors=F)
  
  Mean_counts_relative$sem_max<-Mean_counts_relative$cells_relative_to_ref + (Mean_counts_relative$sd/Mean_counts_relative$samples)
  Mean_counts_relative$sem_min<-Mean_counts_relative$cells_relative_to_ref - (Mean_counts_relative$sd/Mean_counts_relative$samples)
 
  
  cat("Mean_counts_relative_0\n")
  str(Mean_counts_relative)
  cat("\n")
  
  
  ### graph -----
  
  vector_fill<-c(brewer.pal(9, "Greens")[c(6)],brewer.pal(9, "Reds")[c(6)],brewer.pal(9, "Purples")[c(6)],brewer.pal(9, "Blues")[c(5)])
  
  graph_Mean_counts_relative<-ggplot(data=Mean_counts_relative,
                                     aes(x=time_point,
                                         y=cells_relative_to_ref))+          
    geom_line(aes(group=Genotype, 
                  color=Genotype), 
              size=2, alpha=0.8)+
    scale_color_manual(values=vector_fill, drop=F)+      
    geom_point(aes(fill=Genotype, color=Genotype),
               size=3, shape=21, stroke=1, alpha=0.8)+
    geom_errorbar(aes(ymin=sem_min, ymax=sem_max,
                      group=Genotype, 
                      color=Genotype), size=2)+
    scale_color_manual(values=vector_fill, drop=F)+
    scale_fill_manual(values=vector_fill, drop=F)+
    scale_x_discrete(name="time_point", expand = c(0.1, 0.1))+     
    scale_y_continuous(name=paste("Fold increase over", reference,"cells", sep=" "), expand = c(0.1, 0.1))+
    theme_classic()+
    theme(axis.title=element_blank(),
          axis.title.y=element_text(size=12,color="black", family="sans"),
          axis.title.x=element_text(size=12,color="black", family="sans"),
          axis.text.y=element_text(size=8,color="black", family="sans"),
          axis.text.x=element_text(size=8,color="black", family="sans"))+
    theme(legend.title = element_text(size=12),
          legend.text = element_text(size=8),
          legend.key.size = unit(0.5, 'cm'), #change legend key size
          legend.key.height = unit(0.5, 'cm'), #change legend key height
          legend.key.width = unit(0.5, 'cm'), #change legend key width
          legend.position="bottom")+
    guides(fill=guide_legend(nrow=2,byrow=TRUE))
  
  
  
  
  setwd(out)
  
  #svgname<-paste("Heatmap_",paste(selected_annotations, collapse="_"),"_",paste(selected_contrasts, collapse="_"),".svg",sep='')
  
  svgname<-paste("Line_plot_",paste(selected_clone_lines, collapse = "_"),".svg",sep='')
  
  ggsave(svgname,plot=graph_Mean_counts_relative, device ='svg')
  
  write.table(Mean_counts_relative, file=paste("Table_for_line_plots_reference_",reference,".tsv",sep=''), quote=F, row.names = F)
  
  
  
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
    make_option(c("--SeuratObject"), type="character", default=NULL, 
                metavar="type", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--selected_clone_lines"), type="character", default=NULL, 
                metavar="type", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--reference"), type="character", default=NULL, 
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
  
  line_plot_function(opt)

  
  
  
}


###########################################################################

system.time( main() )