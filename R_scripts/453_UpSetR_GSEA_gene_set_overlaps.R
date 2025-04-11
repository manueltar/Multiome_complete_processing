
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
suppressMessages(library("cowplot", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("ggupset", lib.loc = "/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("RColorBrewer", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("splitstackshape", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("org.Hs.eg.db", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))


opt = NULL

options(warn = 1)

Convert_up <- function(TB){
  
  # from https://karobben.github.io/2021/03/31/R/ggupset/ and https://github.com/const-ae/ggupset
  
  TB_tmp = TB
  
  cat("TB_tmp_initial\n")
  cat(str(TB_tmp))
  cat("\n")
  
  for(i in colnames(TB)){
    # cat("-------------------------->i\t")
    # cat(sprintf(as.character(i)))
    # cat("\n")
    
    TB_tmp[i] = i
    
    # cat("TB_tmp_0\n")
    # cat(str(TB_tmp))
    # cat("\n")
    
  }
  TB_tmp[TB == 0] = ""
  
  # cat("TB_tmp_1\n")
  # cat(str(TB_tmp))
  # cat("\n")
  
  TB_t <- data.frame(t(TB_tmp), stringsAsFactors = F)
  
  # cat("TB_t_0\n")
  # cat(str(TB_t))
  # cat("\n")
 
  TB_tmp$upset = as.list(TB_t)
  
  # cat("TB_tmp_FIN\n")
  # cat(str(TB_tmp))
  # cat("\n")
  
  return(TB_tmp)
}

upsetr_function = function(option_list)
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
  
  #### READ and transform indir ----
  
  indir = opt$indir
  
  cat("indir_\n")
  cat(sprintf(as.character(indir)))
  cat("\n")
  
  path_graphs<-paste(out,'graphs/',sep='')
  
  if(file.exists(path_graphs)){
    
  }else{
    
    dir.create(path_graphs)
    
  }#file.exists(path_graphs)
  
  #####################################################
  
  setwd(indir)
  
  files_df<-list.files()
  
  cat("files_df_0\n")
  cat(str(files_df))
  cat("\n")
  
  
  files_df<-list.files()
  
  cat("files_df_0\n")
  cat(str(files_df))
  cat("\n")
  
  files_df_sel<-files_df[grep("Hs\\.entrez_selected\\.rds$",files_df)]
  
  cat("files_df_sel_0\n")
  cat(str(files_df_sel))
  cat("\n")
  
  
  df_files <- data.frame(matrix(NA,
                                nrow = length(files_df_sel),
                                ncol = 2))
  colnames(df_files)[which(colnames(df_files) == 'X1')]<-'collection'
  colnames(df_files)[which(colnames(df_files) == 'X2')]<-'path'
  
  cat("df_files_0\n")
  cat(str(df_files))
  cat("\n")
  
  
  df_files$path<-paste(indir,files_df_sel,sep='')
  df_files$collection<-gsub("_Hs\\.entrez_selected\\.rds$","",files_df_sel)
 
  
  
  cat("df_files_1\n")
  cat(str(df_files))
  cat("\n")
  
 
  ##### Reading loop -----
  
  collection_array<-df_files$collection
  
  cat("collection_array_0\n")
  cat(str(collection_array))
  cat("\n")
  
  
  DEBUG<-0
  
  list_results<-list()
  
  
  for(i in 1:length(collection_array))
  {
    collection_array_sel<-collection_array[i]
    
    cat("------------------------------------------------>\t")
    cat(sprintf(as.character(i)))
    cat("\t")
    cat(sprintf(as.character(collection_array_sel)))
    cat("\n")
    
    df_files_sel<-df_files[which(df_files$collection == collection_array_sel),]
    

    if(DEBUG == 1)
    {
      cat("df_files_sel_0\n")
      cat(str(df_files_sel))
      cat("\n")
      
    }
    
    
    
    if(dim(df_files_sel)[1] >0)
    {
      sel_file<-df_files_sel$path
      
      # cat("--->\t")
      # cat(sprintf(as.character(sel_file)))
      # cat("\n")
      
      SIZE_gate<-file.info(sel_file)$size
      
      if(DEBUG == 1)
      {
        cat("SIZE_gate_0\n")
        cat(str(SIZE_gate))
        cat("\n")
      }
      
      if(SIZE_gate> 0)
      {
        LINE_gate<-length(readLines(sel_file))
        
        if(DEBUG == 1)
        {
          cat("LINE_gate_0\n")
          cat(str(LINE_gate))
          cat("\n")
        }
        
        if(LINE_gate> 0)
        {
          genesets_df<-readRDS(file=sel_file)

          if(DEBUG == 1)
          {
            cat("genesets_df_0\n")
            cat(str(genesets_df))
            cat("\n")
          }
          
          genesets_df$Symbol<-mapIds(org.Hs.eg.db, keys=genesets_df$gene, keytype="ENTREZID",
                                               column="SYMBOL")
          
          if(DEBUG == 1)
          {
            cat("genesets_df_1\n")
            cat(str(genesets_df))
            cat("\n")
          }
          
          ### subgraph -----
          
          subgraph<-genesets_df
          
          subgraph$term<-factor(subgraph$term)
          
          cat("subgraph_0\n")
          cat(str(subgraph))
          cat("\n")
         
          
          
          
          #### subgraph ----
          
          subgraph.dt<-data.table(subgraph,key=c('term'))
          
          cat("subgraph.dt\n")
          cat(str(subgraph.dt))
          cat("\n")
          
          
          subgraph_CLASS_summarised<-as.data.frame(subgraph.dt[,.(instances=.N), by=key(subgraph.dt)],stringsAsFactors=F)
          
          
          
          cat("subgraph_CLASS_summarised\n")
          cat(str(subgraph_CLASS_summarised))
          cat("\n")
          cat(sprintf(as.character(names(summary(subgraph_CLASS_summarised$term)))))
          cat("\n")
          cat(sprintf(as.character(summary(subgraph_CLASS_summarised$term))))
          cat("\n")
          cat(sprintf(as.character(names(summary(subgraph_CLASS_summarised$instances)))))
          cat("\n")
          cat(sprintf(as.character(summary(subgraph_CLASS_summarised$instances))))
          cat("\n")
          
          
          step<-round(max(subgraph_CLASS_summarised$instances)/10,0)
          
          # step<-5
          
          cat("--------->\t")
          cat(sprintf(as.character(step)))
          cat("\n")
          
          if(step == 0){
            
            step=1
          }
          breaks.x<-rev(sort(unique(c(max(subgraph_CLASS_summarised$instances),seq(0,max(subgraph_CLASS_summarised$instances), by=step)))))
          labels.x<-as.character(round(breaks.x,0))
          
          
          cat(sprintf(as.character(breaks.x)))
          cat("\n")
          
          
          subgraph<-ggplot(data=subgraph_CLASS_summarised,
                           aes(y=term,
                               x=instances)) +
            geom_bar(stat="identity",colour='black', fill = "black")+
            theme_bw()+
            theme_classic()+
            theme(axis.title.y=element_blank(),
                  axis.text.y=element_text(angle=0,size=3, color="black", family="sans"),
                  axis.text.x=element_text(angle=0,size=6, color="black", family="sans"),
                  axis.line.x = element_line(size = 0.2),
                  axis.line.y = element_line(size = 0.2),
                  axis.ticks.x = element_line(size = 0.2),
                  axis.ticks.y = element_line(size = 0.2))+
            scale_x_reverse(name="# genes",breaks=breaks.x,labels=labels.x,
                            limits=c(breaks.x[1],breaks.x[length(breaks.x)]))+
            scale_y_discrete(position="right",name=NULL, drop=T)+
            theme(legend.position="hidden")+
            ggeasy::easy_center_title()
          
          cat("subgraph genes DONE\n")
          
          
          #### UpsetR ----
          
          DEBUG<-1
          
          genesets_df$Presence<-1
          
          genesets_df_wide<-unique(as.data.frame(pivot_wider(genesets_df,
                                                                             id_cols=Symbol,
                                                                             names_from=term,
                                                                             values_from=Presence), stringsAsFactors=F))
          
          
          if(DEBUG == 1)
          {
            cat("genesets_df_wide_0\n")
            cat(str(genesets_df_wide))
            cat("\n")
            cat(str(unique(genesets_df_wide$term)))
            cat("\n")
          }
          
          genesets_df_wide[is.na(genesets_df_wide)]<-0
          
          if(DEBUG == 1)
          {
            cat("genesets_df_wide_1\n")
            cat(str(genesets_df_wide))
            cat("\n")
            cat(str(unique(genesets_df_wide$term)))
            cat("\n")
          }
          
          
          
          
          indx.dep<-c(which(colnames(genesets_df_wide) == 'Symbol'))
          
          TMP <- Convert_up(genesets_df_wide[,-indx.dep]) # remove the Symbol & Symbol
          
          if(DEBUG == 1)
          {
            cat("TMP_0\n")
            cat(str(TMP))
            cat("\n")
            
          }
          
          TMP$Symbol<-genesets_df_wide$Symbol

          if(DEBUG == 1)
          {
            cat("TMP_1\n")
            cat(str(TMP))
            cat("\n")
            
          }
          
         
          
          
          
        
          
          #aes(label=after_stat(count)
          
          UpsetR_plot<-ggplot(data=TMP,
                              aes(x = upset)) +
            geom_bar(color='black', fill="steelblue") +
            scale_x_upset(order_by = "freq")+
            theme_combmatrix(
              combmatrix.label.make_space = TRUE,
              combmatrix.label.width = NULL,
              combmatrix.label.height = NULL,
              combmatrix.label.extra_spacing = 3,
              combmatrix.label.total_extra_spacing = unit(10, "pt"),
              combmatrix.label.text = NULL,
              combmatrix.panel.margin = unit(c(1.5, 1.5), "pt"),
              combmatrix.panel.striped_background = TRUE,
              combmatrix.panel.striped_background.color.one = "white",
              combmatrix.panel.striped_background.color.two = "#F7F7F7",
              combmatrix.panel.point.size = 1,
              combmatrix.panel.line.size = 1,
              combmatrix.panel.point.color.fill = "black",
              combmatrix.panel.point.color.empty = "#E0E0E0")
          
          UpsetR_plot<-UpsetR_plot+
            ggtitle(paste("Total gene sets compared =",length(unique(genesets_df$term)),sep=' '))+
            theme_classic()+
            theme(plot.title=element_text(size=8, color="black", family="sans"),
                  axis.title.y=element_text(size=8, color="black", family="sans"),
                  axis.title.x=element_blank(),
                  axis.text.y=element_text(size=8, color="black", family="sans"),
                  axis.text.x=element_blank(),
                  axis.line.x = element_line(size = 0.2),
                  axis.line.y = element_line(size = 0.2),
                  axis.ticks.x = element_line(size = 0.2),
                  axis.ticks.y = element_line(size = 0.2))+
            theme(legend.key.size = unit(0.25, 'cm'), #change legend key size
                  legend.key.height = unit(0.25, 'cm'), #change legend key height
                  legend.key.width = unit(0.25, 'cm'), #change legend key width
                  legend.title = element_text(size=8, family="sans"), #change legend title font size
                  legend.text = element_text(size=6, family="sans"),
                  legend.position="bottom")+ #change legend text font size
            ggeasy::easy_center_title()
          
          
          ##### FINAL arrangement ----
          
          subgraph_FINAL<-plot_grid(NULL,subgraph,
                                    nrow = 2,
                                    ncol=1,
                                    rel_heights = c(0.715, 0.285))
          
          graph_FINAL<-plot_grid(subgraph_FINAL,UpsetR_plot,
                                 nrow = 1,
                                 ncol=2,
                                 rel_widths=c(0.25,0.75))
          
          
          
          
          setwd(path_graphs)
          
          svgname<-paste("UpsetR_","gene_sets_",collection_array_sel,".svg",sep='')
          makesvg = TRUE
          
          if (makesvg == TRUE)
          {
            ggsave(svgname, plot= graph_FINAL,
                   device="svg",
                   height=13, width=13)
          }
          
          
          
        }#LINE_gate> 0
      }# SIZE_gate> 0
    }# dim(df_files_sel)[1] >0
    
  }# i in 1:length(collection_array)
  
 
 
  
  
 
  
  
  
  
 
  
 
}


upsetr_function_restricted = function(option_list)
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
  
  #### READ and transform selected_terms ----
  
  selected_terms = unlist(strsplit(opt$selected_terms, split=","))
  
  cat("selected_terms_\n")
  cat(sprintf(as.character(selected_terms)))
  cat("\n")
  
  
  
  #### READ and transform indir ----
  
  indir = opt$indir
  
  cat("indir_\n")
  cat(sprintf(as.character(indir)))
  cat("\n")
  
  path_graphs<-paste(out,'graphs/',sep='')
  
  if(file.exists(path_graphs)){
    
    unlink(path_graphs)
    dir.create(path_graphs)
    
    
  }else{
    
    dir.create(path_graphs)
    
  }#file.exists(path_graphs)
  
  #####################################################
  
  setwd(indir)
  
  files_df<-list.files()
  
  cat("files_df_0\n")
  cat(str(files_df))
  cat("\n")
  
  
  files_df<-list.files()
  
  cat("files_df_0\n")
  cat(str(files_df))
  cat("\n")
  
  files_df_sel<-files_df[grep("Hs\\.entrez_selected\\.rds$",files_df)]
  
  cat("files_df_sel_0\n")
  cat(str(files_df_sel))
  cat("\n")
  
  
  df_files <- data.frame(matrix(NA,
                                nrow = length(files_df_sel),
                                ncol = 2))
  colnames(df_files)[which(colnames(df_files) == 'X1')]<-'collection'
  colnames(df_files)[which(colnames(df_files) == 'X2')]<-'path'
  
  cat("df_files_0\n")
  cat(str(df_files))
  cat("\n")
  
  
  df_files$path<-paste(indir,files_df_sel,sep='')
  df_files$collection<-gsub("_Hs\\.entrez_selected\\.rds$","",files_df_sel)
  
  
  
  cat("df_files_1\n")
  cat(str(df_files))
  cat("\n")
  
  
  ##### Reading loop -----
  
  collection_array<-df_files$collection
  
  cat("collection_array_0\n")
  cat(str(collection_array))
  cat("\n")
  
  
  DEBUG<-1
  
  list_results<-list()
  
  
  for(i in 1:length(collection_array))
  {
    collection_array_sel<-collection_array[i]
    
    cat("------------------------------------------------>\t")
    cat(sprintf(as.character(i)))
    cat("\t")
    cat(sprintf(as.character(collection_array_sel)))
    cat("\n")
    
    df_files_sel<-df_files[which(df_files$collection == collection_array_sel),]
    
    
    if(DEBUG == 1)
    {
      cat("df_files_sel_0\n")
      cat(str(df_files_sel))
      cat("\n")
      
    }
    
    
    
    if(dim(df_files_sel)[1] >0)
    {
      sel_file<-df_files_sel$path
      
      # cat("--->\t")
      # cat(sprintf(as.character(sel_file)))
      # cat("\n")
      
      SIZE_gate<-file.info(sel_file)$size
      
      if(DEBUG == 1)
      {
        cat("SIZE_gate_0\n")
        cat(str(SIZE_gate))
        cat("\n")
      }
      
      if(SIZE_gate> 0)
      {
        LINE_gate<-length(readLines(sel_file))
        
        if(DEBUG == 1)
        {
          cat("LINE_gate_0\n")
          cat(str(LINE_gate))
          cat("\n")
        }
        
        if(LINE_gate> 0)
        {
          genesets_df<-readRDS(file=sel_file)
          
          if(DEBUG == 1)
          {
            cat("genesets_df_0\n")
            cat(str(genesets_df))
            cat("\n")
          }
          
          genesets_df$Symbol<-mapIds(org.Hs.eg.db, keys=genesets_df$gene, keytype="ENTREZID",
                                     column="SYMBOL")
          
          if(DEBUG == 1)
          {
            cat("genesets_df_1\n")
            cat(str(genesets_df))
            cat("\n")
          }
          
          matches <- grep(paste(selected_terms,collapse="|"),genesets_df$term)
          
          toMatch<-tolower(selected_terms)
          
          
          matches_lc <- grep(paste(toMatch,collapse="|"),genesets_df$term)
          
          total_matches<-unique(c(matches,matches_lc))
          
          
          if(length(total_matches) >0){
            
            genesets_df_sel<-genesets_df[total_matches,]
            
            cat("genesets_df_sel_0\n")
            cat(str(genesets_df_sel))
            cat("\n")
            
            ### subgraph -----
            
            
            subgraph<-genesets_df_sel
            
            subgraph$term<-factor(subgraph$term)
            
            if(DEBUG == 1)
            {
              cat("subgraph_0\n")
              cat(str(subgraph))
              cat("\n")
            }
            


            

            subgraph.dt<-data.table(subgraph,key=c('term'))
            
            if(DEBUG == 1)
            {
              cat("subgraph.dt\n")
              cat(str(subgraph.dt))
              cat("\n")
            }
            
            
            subgraph_CLASS_summarised<-as.data.frame(subgraph.dt[,.(instances=.N), by=key(subgraph.dt)],stringsAsFactors=F)
            
            
            if(DEBUG == 1)
            {
              cat("subgraph_CLASS_summarised\n")
              cat(str(subgraph_CLASS_summarised))
              cat("\n")
              cat(sprintf(as.character(names(summary(subgraph_CLASS_summarised$term)))))
              cat("\n")
              cat(sprintf(as.character(summary(subgraph_CLASS_summarised$term))))
              cat("\n")
              cat(sprintf(as.character(names(summary(subgraph_CLASS_summarised$instances)))))
              cat("\n")
              cat(sprintf(as.character(summary(subgraph_CLASS_summarised$instances))))
              cat("\n")
            }
            
            
            step<-round(max(subgraph_CLASS_summarised$instances)/10,0)
            
            # step<-5
            
            cat("--------->\t")
            cat(sprintf(as.character(step)))
            cat("\n")
            
            if(step == 0){
              
              step=1
            }
            breaks.x<-rev(sort(unique(c(max(subgraph_CLASS_summarised$instances),seq(0,max(subgraph_CLASS_summarised$instances), by=step)))))
            labels.x<-as.character(round(breaks.x,0))
            
            if(DEBUG == 1)
            {
              cat(sprintf(as.character(breaks.x)))
              cat("\n")
            }
            
            
            subgraph<-ggplot(data=subgraph_CLASS_summarised,
                             aes(y=term,
                                 x=instances)) +
              geom_bar(stat="identity",colour='black', fill = "black")+
              theme_bw()+
              theme_classic()+
              theme(axis.title.y=element_blank(),
                    axis.text.y=element_text(angle=0,size=3, color="black", family="sans"),
                    axis.text.x=element_text(angle=0,size=6, color="black", family="sans"),
                    axis.line.x = element_line(size = 0.2),
                    axis.line.y = element_line(size = 0.2),
                    axis.ticks.x = element_line(size = 0.2),
                    axis.ticks.y = element_line(size = 0.2))+
              scale_x_reverse(name="# genes",breaks=breaks.x,labels=labels.x,
                              limits=c(breaks.x[1],breaks.x[length(breaks.x)]))+
              scale_y_discrete(position="right",name=NULL, drop=T)+
              theme(legend.position="hidden")+
              ggeasy::easy_center_title()
            
            cat("subgraph genes DONE\n")
            
            
            #### UpsetR ----
            

            genesets_df_sel$Presence<-1
            
            genesets_df_sel_wide<-unique(as.data.frame(pivot_wider(genesets_df_sel,
                                                               id_cols=Symbol,
                                                               names_from=term,
                                                               values_from=Presence), stringsAsFactors=F))
            
            
            if(DEBUG == 1)
            {
              cat("genesets_df_sel_wide_0\n")
              cat(str(genesets_df_sel_wide))
              cat("\n")
              cat(str(unique(genesets_df_sel_wide$term)))
              cat("\n")
            }
            
            genesets_df_sel_wide[is.na(genesets_df_sel_wide)]<-0
            
            if(DEBUG == 1)
            {
              cat("genesets_df_sel_wide_1\n")
              cat(str(genesets_df_sel_wide))
              cat("\n")
              cat(str(unique(genesets_df_sel_wide$term)))
              cat("\n")
            }
            
            
            
            
            indx.dep<-c(which(colnames(genesets_df_sel_wide) == 'Symbol'))
            
            TMP <- Convert_up(genesets_df_sel_wide[,-indx.dep]) # remove the Symbol & Symbol
            
            if(DEBUG == 1)
            {
              cat("TMP_0\n")
              cat(str(TMP))
              cat("\n")
              
            }
            
            TMP$Symbol<-genesets_df_sel_wide$Symbol
            
            if(DEBUG == 1)
            {
              cat("TMP_1\n")
              cat(str(TMP))
              cat("\n")
              
            }
            
            
            
            
            
            
            
            #aes(label=after_stat(count)
            
            UpsetR_plot<-ggplot(data=TMP,
                                aes(x = upset)) +
              geom_bar(color='black', fill="steelblue") +
              scale_x_upset(order_by = "freq")+
              theme_combmatrix(
                combmatrix.label.make_space = TRUE,
                combmatrix.label.width = NULL,
                combmatrix.label.height = NULL,
                combmatrix.label.extra_spacing = 3,
                combmatrix.label.total_extra_spacing = unit(10, "pt"),
                combmatrix.label.text = NULL,
                combmatrix.panel.margin = unit(c(1.5, 1.5), "pt"),
                combmatrix.panel.striped_background = TRUE,
                combmatrix.panel.striped_background.color.one = "white",
                combmatrix.panel.striped_background.color.two = "#F7F7F7",
                combmatrix.panel.point.size = 1,
                combmatrix.panel.line.size = 1,
                combmatrix.panel.point.color.fill = "black",
                combmatrix.panel.point.color.empty = "#E0E0E0")
            
            UpsetR_plot<-UpsetR_plot+
              ggtitle(paste("Total gene sets compared =",length(unique(genesets_df_sel$term)),sep=' '))+
              theme_classic()+
              theme(plot.title=element_text(size=8, color="black", family="sans"),
                    axis.title.y=element_text(size=8, color="black", family="sans"),
                    axis.title.x=element_blank(),
                    axis.text.y=element_text(size=8, color="black", family="sans"),
                    axis.text.x=element_blank(),
                    axis.line.x = element_line(size = 0.2),
                    axis.line.y = element_line(size = 0.2),
                    axis.ticks.x = element_line(size = 0.2),
                    axis.ticks.y = element_line(size = 0.2))+
              theme(legend.key.size = unit(0.25, 'cm'), #change legend key size
                    legend.key.height = unit(0.25, 'cm'), #change legend key height
                    legend.key.width = unit(0.25, 'cm'), #change legend key width
                    legend.title = element_text(size=8, family="sans"), #change legend title font size
                    legend.text = element_text(size=6, family="sans"),
                    legend.position="bottom")+ #change legend text font size
              ggeasy::easy_center_title()
            
            
            ##### FINAL arrangement ----
            
            subgraph_FINAL<-plot_grid(NULL,subgraph,
                                      nrow = 2,
                                      ncol=1,
                                      rel_heights = c(0.715, 0.285))
            
            graph_FINAL<-plot_grid(subgraph_FINAL,UpsetR_plot,
                                   nrow = 1,
                                   ncol=2,
                                   rel_widths=c(0.25,0.75))
            
            
            
            
            setwd(path_graphs)
            
            svgname<-paste("Selected_UpsetR_","gene_sets_",collection_array_sel,".svg",sep='')
            makesvg = TRUE
            
            if (makesvg == TRUE)
            {
              ggsave(svgname, plot= graph_FINAL,
                     device="svg",
                     height=13, width=13)
            }
            
            
            
          }#length(total_matches) >0
        }#LINE_gate> 0
      }# SIZE_gate> 0
    }# dim(df_files_sel)[1] >0
    
  }# i in 1:length(collection_array)
  
  
  
  
  
  
  
  
  
  
  
  
  
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
    make_option(c("--indir"), type="character", default=NULL, 
                metavar="type", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--selected_terms"), type="character", default=NULL, 
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
  
  # upsetr_function(opt)
  upsetr_function_restricted(opt)

  
}


###########################################################################

system.time( main() )