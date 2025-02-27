
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


opt = NULL

options(warn = 1)




AP_function = function(option_list)
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
  
  #### READ and transform out ----
  
  exemption_terms = unlist(strsplit(opt$exemption_terms, split=","))
  
  cat("exemption_terms_0\n")
  cat(sprintf(as.character(exemption_terms)))
  cat("\n")
  
 
  #### READ and transform minGS_size ----
  
  minGS_size = opt$minGS_size
  
  cat("minGS_size_\n")
  cat(sprintf(as.character(minGS_size)))
  cat("\n")
  
  #### READ and transform maxGS_size ----
  
  maxGS_size = opt$maxGS_size
  
  cat("maxGS_size_\n")
  cat(sprintf(as.character(maxGS_size)))
  cat("\n")
  
  #### READ and transform pval_threshold ----
  
  pval_threshold = opt$pval_threshold
  
  cat("pval_threshold_\n")
  cat(sprintf(as.character(pval_threshold)))
  cat("\n")
  
  #### READ and transform log2FC_threshold ----
  
  log2FC_threshold = opt$log2FC_threshold
  
  cat("log2FC_threshold_\n")
  cat(sprintf(as.character(log2FC_threshold)))
  cat("\n")
  
  #### READ DE_results ----
  
 
  DE_results<-as.data.frame(fread(file=opt$DE_results,sep="\t", header=TRUE), stringsAsFactors=F)
  
  cat("DE_results_0\n")
  cat(str(DE_results))
  cat("\n")
  cat(sprintf(as.character(names(summary(as.factor(DE_results$contrast))))))
  cat("\n")
  cat(sprintf(as.character(summary(as.factor(DE_results$contrast)))))
  cat("\n")
  
  
 
  
  DE_results$ENTREZID <- mapIds(org.Hs.eg.db, keys=DE_results$gene, keytype="SYMBOL",
                                 column="ENTREZID")
  
  # cat("DE_results_1\n")
  # cat(str(DE_results))
  # cat("\n")
  
  
  DE_results_NO_NA<-DE_results[!is.na(DE_results$ENTREZID) &
                                                                !is.na(DE_results$padj),]
  
  
  DEBUG<-0
  
  LIST_TRUE_FINAL<-list()
  
  ALL_tested_TRUE_FINAL<-list()
  
  if(dim(DE_results_NO_NA)[1] >0)
  {
  
    cat("DE_results_NO_NA_0\n")
    cat(str(DE_results_NO_NA))
    cat("\n")
    cat(sprintf(as.character(names(summary(as.factor(DE_results_NO_NA$contrast))))))
    cat("\n")
    cat(sprintf(as.character(summary(as.factor(DE_results_NO_NA$contrast)))))
    cat("\n")
    
    array_contrasts<-unique(DE_results_NO_NA$contrast)
    
    cat("array_contrasts_0\n")
    cat(str(array_contrasts))
    cat("\n")
    
    
    for(iteration_array_contrasts in 1:length(array_contrasts))
    {
      array_contrasts_sel<-array_contrasts[iteration_array_contrasts]
      
      cat("-------------------------------------------------------------------------------------------------------------------------------------------------------->\t")
      cat(sprintf(as.character(array_contrasts_sel)))
      cat("\n")
      
      DE_results_NO_NA_sel<-DE_results_NO_NA[which(DE_results_NO_NA$contrast == array_contrasts_sel),]
      
      
      if(DEBUG ==1){
        
        cat("DE_results_NO_NA_sel_0\n")
        cat(str(DE_results_NO_NA_sel))
        cat("\n")
        cat(sprintf(as.character(names(summary(as.factor(DE_results_NO_NA_sel$contrast))))))
        cat("\n")
        cat(sprintf(as.character(summary(as.factor(DE_results_NO_NA_sel$contrast)))))
        cat("\n")
      }
      
      array_cell_types<-unique(DE_results_NO_NA_sel$cell_type)
      
      if(DEBUG ==1){
      
        
        cat("array_cell_types_0\n")
        cat(str(array_cell_types))
        cat("\n")
      }
      
      List_cell_type<-list()
      
      for(iteration_2 in 1:length(array_cell_types)){
        
        array_cell_types_sel<-array_cell_types[iteration_2]
        
        cat("------------------------------------------>\t")
        cat(sprintf(as.character(array_cell_types_sel)))
        cat("\n")
        

        DE_results_tmp<-DE_results_NO_NA_sel[which(DE_results_NO_NA_sel$cell_type == array_cell_types_sel),]
        
        
        if(DEBUG ==1){
          
          cat("DE_results_tmp_0\n")
          cat(str(DE_results_tmp))
          cat("\n")
          cat(sprintf(as.character(names(summary(as.factor(DE_results_tmp$contrast))))))
          cat("\n")
          cat(sprintf(as.character(summary(as.factor(DE_results_tmp$contrast)))))
          cat("\n")
        }
        
        
        indx.DUP<-which(duplicated(DE_results_tmp$ENTREZID) == TRUE)
        
        if(DEBUG ==1){
          cat("indx.DUP_0\n")
          cat(str(indx.DUP))
          cat("\n")
        }
        
        if(length(indx.DUP) >0)
        {
          
          DE_results_tmp_NO_DUP<-DE_results_tmp[-indx.DUP,]
          
        }else{
          
          DE_results_tmp_NO_DUP<-DE_results_tmp
          
        }#length(indx.DUP) >0
        
        
        
        if(dim(DE_results_tmp_NO_DUP)[1] >0)
        {
          DE_results_tmp_NO_DUP$AbsLogC<-abs(DE_results_tmp_NO_DUP$log2FoldChange)
          
          if(DEBUG ==1){
            cat("DE_results_tmp_NO_DUP_0\n")
            cat(str(DE_results_tmp_NO_DUP))
            cat("\n")
          }
          
          
          DE_results_tmp_NO_DUP_FC_thresholded<-DE_results_tmp_NO_DUP[which(DE_results_tmp_NO_DUP$AbsLogC >= log2FC_threshold),]
          
          
          if(dim(DE_results_tmp_NO_DUP_FC_thresholded)[1] >0)
          {
            if(DEBUG ==1){
              cat("DE_results_tmp_NO_DUP_FC_thresholded_0\n")
              cat(str(DE_results_tmp_NO_DUP_FC_thresholded))
              cat("\n")
            }
            
            
            check.DE<-DE_results_tmp_NO_DUP_FC_thresholded[which(DE_results_tmp_NO_DUP_FC_thresholded$padj <= pval_threshold),]
            
            
            if(DEBUG ==1){
              cat("check.DE_0\n")
              cat(str(check.DE))
              cat("\n")
            }
            
            if(dim(check.DE) [1] >1)
            {
              
              
              # ###############################################################
              # quit(status = 1)
              
              
              scores<-matrix(DE_results_tmp_NO_DUP_FC_thresholded$padj)
              row.names(scores)<-DE_results_tmp_NO_DUP_FC_thresholded$ENTREZID
              colnames(scores)<-"padj"
              
              if(DEBUG ==1){
                cat("scores_0\n")
                cat(str(scores))
                cat("\n")
              }
              
              #### READ and transform background_genes ----
              
              background_genes = opt$background_genes
              
              HPA<-read.GMT(background_genes)
              
              # cat("HPA_0\n")
              # cat(str(HPA))
              # cat("\n")
              # 
              background<-makeBackground(HPA)
              
              if(DEBUG ==1){
                cat("background_0\n")
                cat(str(background))
                cat("\n")
              }
              
              #### READ and transform path_to_GMT ----
              
              path_to_GMT = opt$path_to_GMT
              
              if(DEBUG ==1){
                cat("path_to_GMT_\n")
                cat(sprintf(as.character(path_to_GMT)))
                cat("\n")
              }
              
              #### READ and transform search_terms ----
              
              search_terms = unlist(strsplit(opt$search_terms, split=","))
              
              if(DEBUG ==1){
                cat("search_terms_\n")
                cat(sprintf(as.character(search_terms)))
                cat("\n")
              }
              
              
              #### list all files in path_to_GMT ----
              
              file_list <- list.files(path=path_to_GMT, include.dirs = FALSE)
              
              if(DEBUG ==1){
                cat("file_list\n")
                cat(str(file_list))
                cat("\n")
              }
              
              
              indexes_sel <- grep("Hs\\.entrez\\.gmt$",file_list)
              
              file_list_sel <- as.data.frame(file_list[indexes_sel], stringsAsFactors=F)
              colnames(file_list_sel)<-"file"
              
              if(DEBUG ==1){
                cat("file_list_sel\n")
                cat(str(file_list_sel))
                cat("\n")
              }
              
              ### Loop to open gmt files ---
              
              # START<-grep('c5.hpo.v2023.2.Hs.entrez.gmt',file_list_sel$file)
              
              START<-1
              
              if(DEBUG ==1){
                cat("START\n")
                cat(str(START))
                cat("\n")
              }
              
              tested<-list()
              
              # DEBUG<-0
              
              setwd(path_to_GMT)
              
              Final_result<-list()
              
              for(i in START:dim(file_list_sel)[1])
              {
                file_sel<-file_list_sel$file[i]
                
                cat("------>\t")
                cat(sprintf(as.character(i)))
                cat("\t")
                cat(sprintf(as.character(file_sel)))
                cat("\n")
                
                
                gmt_sel<-read.GMT(file_sel)
                
                # cat(sprintf(as.character(length(gmt_sel))))
                # cat("\n")
                
                
                # if(gmt_sel == 'c5.hpo.v2023.2.Hs.entrez.gmt')
                # {
                #   DEBUG<-1 
                # }else{
                #   
                #   DEBUG<-0
                # }
                
                # if(DEBUG == 1)
                # {
                #   cat("gmt_sel_0\n")
                #   cat(str(gmt_sel))
                #   cat("\n")
                # }
                
                id_vector<-vector()
                
                for(k in 1:length(gmt_sel))
                {
                  id_list<-unique(gmt_sel[[k]]$id)
                  
                  # if(DEBUG == 1)
                  # {
                  #   cat("id_list_0\n")
                  #   cat(str(id_list))
                  #   cat("\n")
                  #   
                  # }
                  
                  id_vector[k]<-paste(id_list, collapse='__')
                  
                  # if(DEBUG == 1)
                  # {
                  #   cat("id_vector_0\n")
                  #   cat(str(id_vector))
                  #   cat("\n")
                  #   
                  # }
                  
                }#k in 1:length(gmt_sel)
                
                if(DEBUG == 1)
                {
                  cat("id_vector_0\n")
                  cat(str(id_vector))
                  cat("\n")
                }
                
                # toMatch<-c("Chr1p13","cHr1p21")
                
                toMatch<-search_terms
                
                matches <- grep(paste(toMatch,collapse="|"),id_vector)
                
                if(DEBUG == 1)
                {
                  cat("matches_0\n")
                  cat(str(matches))
                  cat("\n")
                }
                
                toMatch<-tolower(search_terms)
                # toMatch<-tolower(c("Chr1p13","cHr1p21"))
                
                matches_lc <- grep(paste(toMatch,collapse="|"),id_vector)
                
                if(DEBUG == 1)
                {
                  cat("matches_lc_0\n")
                  cat(str(matches_lc))
                  cat("\n")
                }
                
                
                total_matches<-unique(c(matches,matches_lc))
                
                gmt_sel_GREP<-gmt_sel[total_matches]
                
                
                
                if(DEBUG == 1)
                {
                  cat("gmt_sel_GREP_0\n")
                  cat(str(gmt_sel_GREP))
                  cat("\n")
                }
                
                if(length(gmt_sel_GREP) >0)
                {
                  
                  Selected_pathways<-data.frame()
                  for(iteration_gmt_sel_GREP in 1:length(gmt_sel_GREP))
                  {
                    reference_gmt_sel<-unique(as.data.frame(gmt_sel_GREP[[iteration_gmt_sel_GREP]], stringsAsFactors=F))
                    
                    if(DEBUG ==1)
                    {
                      cat("reference_gmt_sel_0\n")
                      cat(str(reference_gmt_sel))
                      cat("\n")
                    }
                    
                    if(dim(reference_gmt_sel)[1] >= minGS_size & dim(reference_gmt_sel)[1] <= maxGS_size)
                    {
                      colnames(reference_gmt_sel)[which(colnames(reference_gmt_sel) == 'genes')]<-'ENTREZ'
                      
                      
                      Selected_pathways<-rbind(reference_gmt_sel,Selected_pathways)
                    }else{
                      
                      #Include Dorothea gene sets
                      
                      # gmt_sel_GREP
                      
                      indx<-which(unique(reference_gmt_sel$id)%in%exemption_terms)
                      
                      FLAG_oversized_but_interesting<-length(indx)
                      
                      
                      if(DEBUG ==1)
                      {
                        cat("FLAG_oversized_but_interesting_0\n")
                        cat(str(FLAG_oversized_but_interesting))
                        cat("\n")
                      }
                      
                      if(FLAG_oversized_but_interesting > 0)
                      {
                        colnames(reference_gmt_sel)[which(colnames(reference_gmt_sel) == 'genes')]<-'ENTREZ'
                        Selected_pathways<-rbind(reference_gmt_sel,Selected_pathways)
                        
                        
                        cat("Hello_world_FLAG_oversized_but_interesting\n")
                        cat(sprintf(as.character(unique(reference_gmt_sel$name))))
                        cat("\n")
                        
                      }#FLAG_oversized_but_interesting > 0
                      
                    }# dim(reference_gmt_sel)[1] >= minGS_size & dim(reference_gmt_sel)[1] <= maxGS_size
                  }#iteration_gmt_sel_GREP in 1:length(gmt_sel_GREP)
                  
                  if(DEBUG ==1)
                  {
                    cat("Selected_pathways_0\n")
                    cat(str(Selected_pathways))
                    cat("\n")
                  }
                  
                  ALL_names<-names(gmt_sel_GREP)
                  
                  if(DEBUG ==1)
                  {
                    cat("ALL_names_0\n")
                    cat(str(ALL_names))
                    cat("\n")
                  }
                  
                  Selected_names<-unique(Selected_pathways$id)
                  
                  if(DEBUG ==1)
                  {
                    cat("Selected_names_0\n")
                    cat(str(Selected_names))
                    cat("\n")
                    cat(sprintf(as.character(Selected_names)))
                    cat("\n")
                  }
                  
                  indx.select<-which(ALL_names%in%Selected_names)
                  
                  if(DEBUG ==1)
                  {
                    cat("indx.select_0\n")
                    cat(str(indx.select))
                    cat("\n")
                  }
                  
                  assayed_gmt<-gmt_sel_GREP[indx.select]
                  
                  if(DEBUG ==1)
                  {
                    cat("assayed_gmt\n")
                    cat(str(assayed_gmt))
                    cat("\n")
                  }
                  
                  assayed_name_vector<-NULL
                  for(iteration_ALL_tested in 1:length(assayed_gmt)){
                    
                    assayed_gmt_sel<-assayed_gmt[[iteration_ALL_tested]]
                    
                    # if(DEBUG ==1)
                    # {
                    #   cat("assayed_gmt_sel\n")
                    #   cat(str(assayed_gmt_sel))
                    #   cat("\n")
                    # }
                    
                    assayed_gmt_sel_name<-assayed_gmt_sel$name
                    
                    
                    if(DEBUG ==1)
                    {
                      cat("assayed_gmt_sel_name\n")
                      cat(str(assayed_gmt_sel_name))
                      cat("\n")
                    }
                    
                    assayed_name_vector[iteration_ALL_tested]<-assayed_gmt_sel_name
                    
                  }#iteration_ALL_tested in 1:length(assayed_gmt)
                  
                  if(DEBUG ==1)
                  {
                    cat("assayed_name_vector_0\n")
                    cat(str(assayed_name_vector))
                    cat("\n")
                  }
                  
                  
                  tested[[i]]<-assayed_name_vector
                  
                  ####ActivePathways seurat_cell_type----
                  
                  AP_result<-ActivePathways(scores, assayed_gmt, background= background,
                                            cutoff = pval_threshold,
                                            significant = 0.05)
                  
                  
                  FLAG_null<-sum(is.null(AP_result))
                  
                  if(FLAG_null == 0)
                  {
                    
                    if(DEBUG == 1)
                    {
                      cat("AP_result_0\n")
                      cat(str(AP_result))
                      cat("\n")
                    }
                    
                    AP_result_df<-data.frame()
                    
                    for(iteration_AP_result in 1:dim(AP_result)[1])
                    {
                      AP_result_sel<-AP_result[iteration_AP_result,]
                      
                      if(DEBUG == 1)
                      {
                        cat("AP_result_sel_0\n")
                        cat(str(AP_result_sel))
                        cat("\n")
                      }
                      
                      ENTREZ_vector<-unique(as.character(unlist(AP_result_sel$overlap)))
                      
                      if(DEBUG == 1)
                      {
                        cat("ENTREZ_vector_0\n")
                        cat(str(ENTREZ_vector))
                        cat("\n")
                      }
                      
                      a.df<-as.data.frame(cbind(rep(as.character(AP_result_sel$term.id), length(ENTREZ_vector)),
                                                rep(as.character(AP_result_sel$term.name), length(ENTREZ_vector)),
                                                rep(as.character(AP_result_sel$adjusted.p.val), length(ENTREZ_vector)),
                                                rep(as.character(AP_result_sel$term.size), length(ENTREZ_vector)),
                                                ENTREZ_vector),stringsAsFactors=F)
                      
                      if(DEBUG == 1)
                      {
                        cat("a.df_0\n")
                        cat(str(a.df))
                        cat("\n")
                      }
                      
                      colnames(a.df)<-c("term.id","term.name","adjusted.p.val","term.size","ENTREZ")
                      
                      a.df$adjusted.p.val<-as.numeric(a.df$adjusted.p.val)
                      a.df$term.size<-as.numeric(a.df$term.size)
                      
                      if(DEBUG == 1)
                      {
                        cat("a.df_1\n")
                        cat(str(a.df))
                        cat("\n")
                      }
                      
                      AP_result_df<-unique(rbind(a.df,AP_result_df))
                      
                    }#iteration_AP_result in 1:dim(AP_result)[1]
                    
                    
                    if(dim(AP_result_df)[1] >0)
                    {
                      if(DEBUG == 1)
                      {
                        cat("AP_result_df_0\n")
                        cat(str(AP_result_df))
                        cat("\n")
                      }
                      
                      colnames(AP_result_df)[which(colnames(AP_result_df) == 'term.id')]<-'id'
                      colnames(AP_result_df)[which(colnames(AP_result_df) == 'term.name')]<-'name'
                      colnames(AP_result_df)[which(colnames(AP_result_df) == 'term.size')]<-'size'
                      colnames(AP_result_df)[which(colnames(AP_result_df) == 'overlap')]<-'ENTREZ'
                      
                      
                      
                      AP_result_df<-merge(AP_result_df,
                                          Selected_pathways,
                                          by=c('id','name','ENTREZ'))
                      
                      if(DEBUG == 1)
                      {
                        cat("AP_result_df_1\n")
                        cat(str(AP_result_df))
                        cat("\n")
                      }
                      
                      AP_result_df$collection<-file_sel
                      
                      if(DEBUG == 1)
                      {
                        cat("AP_result_df_2\n")
                        cat(str(AP_result_df))
                        cat("\n")
                      }
                      
                      AP_result_df$ensembl_gene_id<-NA
                      AP_result_df$ensembl_gene_id<-mapIds(org.Hs.eg.db, keys=AP_result_df$ENTREZ, keytype="ENTREZID",
                                                           column="ENSEMBL")
                      
                      
                      
                      multiVals <- function(x) paste(x,collapse=";")
                      
                      AP_result_df$Symbol<-NA
                      
                      AP_result_df$Symbol <- mapIds(org.Hs.eg.db, keys=AP_result_df$ENTREZ, keytype="ENTREZID",
                                                    column="SYMBOL", multiVals=multiVals)
                      
                      if(DEBUG ==1)
                      {
                        cat("AP_result_df_3\n")
                        cat(str(AP_result_df))
                        cat("\n")
                      }
                      
                      AP_result_df.dt<-data.table(AP_result_df, key=c('collection','id','name','adjusted.p.val','size'))
                      
                      
                      AP_result_df_collapsed<-as.data.frame(AP_result_df.dt[,.(n_genes_in_overlap=length(ensembl_gene_id),
                                                                               overlap_symbol=paste(Symbol, collapse='|'),
                                                                               overlap_ensembl_gene_id=paste(ensembl_gene_id, collapse='|')), by=key(AP_result_df.dt)])
                      
                      if(DEBUG ==1)
                      {
                        cat("AP_result_df_collapsed_0\n")
                        cat(str(AP_result_df_collapsed))
                        cat("\n")
                      }
                      
                      # check<-AP_result_df_collapsed[which(AP_result_df_collapsed$id == 'HP_ABNORMAL_PLATELET_VOLUME'),]
                      # 
                      # if(DEBUG ==1)
                      # {
                      #   cat("check_0\n")
                      #   cat(str(check))
                      #   cat("\n")
                      #   cat(sprintf(as.character(check$overlap_symbol)))
                      #   cat("\n")
                      # }
                      
                      Final_result[[i]]<-AP_result_df_collapsed
                      
                    }#dim(AP_result_df)[1] >0
                  }#FLAG_null == 0
                }#length(gmt_sel_GREP) >1
              }# i in 1:dim(file_list_sel)[1]
              
              if(length(tested) >0){
                
                ALL_tested = unique(unlist(tested))
                
                cat("ALL_tested_0\n")
                cat(str(ALL_tested))
                cat("\n")
                
                
                ALL_tested_sel.df<-as.data.frame(ALL_tested)
                
                colnames(ALL_tested_sel.df)<-"Tested_Pathways"
                
                cat("ALL_tested_sel.df_0\n")
                cat(str(ALL_tested_sel.df))
                cat("\n")
                
                
                ALL_tested_TRUE_FINAL[[iteration_array_contrasts]]<-ALL_tested_sel.df
                
                
                
              }# length(tested) >0
              
              if(length(Final_result)>0)
              {
                FINAL_df = unique(as.data.frame(data.table::rbindlist(Final_result, fill = T)))
                
                cat("FINAL_df_0\n")
                cat(str(FINAL_df))
                cat("\n")
                
                FINAL_df$Minus_logpval<-round(-1*log10(FINAL_df$adjusted.p.val),2)
                
                cat("FINAL_df_1\n")
                cat(str(FINAL_df))
                cat("\n")
                
                FINAL_df$Significance<-NA
                
                FINAL_df$Significance[which(FINAL_df$Minus_logpval >= 1.3)]<-'YES'
                FINAL_df$Significance[which(FINAL_df$Minus_logpval < 1.3)]<-'NO'
                
                FINAL_df$Significance<-factor(FINAL_df$Significance,
                                              levels=c('NO','YES'))
                
                cat("FINAL_df_2\n")
                cat(str(FINAL_df))
                cat("\n")
                
                FINAL_df.dt<-data.table(FINAL_df, key='id')
                
                
                
                FINAL_df_MAX<-as.data.frame(FINAL_df.dt[,.SD[which.max(Minus_logpval)], by=key(FINAL_df.dt)], stringsAsFactors=F)
                
                FINAL_df_MAX$cell_type<-array_cell_types_sel

                cat("FINAL_df_MAX_0\n")
                cat(str(FINAL_df_MAX))
                cat("\n")
                cat(str(unique(FINAL_df_MAX$id)))
                cat("\n")
                
                # ALL_tested = unique(as.data.frame(data.table::rbindlist(tested, fill = T)))
                
                # cat("ALL_tested_0\n")
                # cat(str(ALL_tested))
                # cat("\n")
                
                List_cell_type[[iteration_2]]<-FINAL_df_MAX
                
              }#length(Final_result)>0
              
              
              
              
            }#dim(check.DE) [1] >0
          }#dim(DE_results_tmp_NO_DUP_FC_thresholded)[1] >0
        }#dim(DE_results_tmp_NO_DUP)[1] >0
        
      }#iteration_2 in 1:length(array_cell_types)
      
      
      if(length(List_cell_type)>0)
      {
        intermediate = unique(as.data.frame(data.table::rbindlist(List_cell_type, fill = T)))
        
        cat("intermediate_0\n")
        cat(str(intermediate))
        cat("\n")
        
        intermediate$contrast<-array_contrasts_sel
        
        
        LIST_TRUE_FINAL[[iteration_array_contrasts]]<-intermediate
        
      }#length(List_cell_type)>0
      
    }#iteration_array_contrasts in 1:length(array_contrasts)
    
    if(length(ALL_tested_TRUE_FINAL)>0)
    {
      BUCEPHALUS = unique(as.data.frame(data.table::rbindlist(ALL_tested_TRUE_FINAL, fill = T)))
      
      cat("BUCEPHALUS_0\n")
      cat(str(BUCEPHALUS))
      cat("\n")
     
      
      path_ORA<-paste(out,'ORA','/',sep='')
      
      if (file.exists(path_ORA)){
        
        unlink(path_ORA, recursive =T)
        dir.create(file.path(path_ORA))
        
        
        
      }else{
        
        dir.create(file.path(path_ORA))
        
      }#path_ORA
      
      setwd(path_ORA)
      
      write.table(BUCEPHALUS,file=paste("ALL_tested_ORA_MSigDB",".tsv", sep=''),sep="\t", quote=F, row.names = F)
      
    }#length(ALL_tested_TRUE_FINAL)>0
    
    if(length(LIST_TRUE_FINAL)>0)
    {
      ESTATIRA = unique(as.data.frame(data.table::rbindlist(LIST_TRUE_FINAL, fill = T)))
      
      cat("ESTATIRA_0\n")
      cat(str(ESTATIRA))
      cat("\n")
      cat(sprintf(as.character(names(summary(as.factor(ESTATIRA$contrast))))))
      cat("\n")
      cat(sprintf(as.character(summary(as.factor(ESTATIRA$contrast)))))
      cat("\n")
      
      path_ORA<-paste(out,'ORA','/',sep='')
      
      if (file.exists(path_ORA)){
        
      }else{
        
        dir.create(file.path(path_ORA))
        
      }#path_ORA
      
      setwd(path_ORA)
      
      write.table(ESTATIRA,file=paste("ORA_MSigDB_results",".tsv", sep=''),sep="\t", quote=F, row.names = F)
      
    }#length(LIST_TRUE_FINAL)>0
  }# dim(DE_results_NO_NA)[1] >0
  
 
  
 
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
    make_option(c("--path_to_GMT"), type="character", default=NULL, 
                metavar="type", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--exemption_terms"), type="character", default=NULL, 
                metavar="type", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--search_terms"), type="character", default=NULL, 
                metavar="type", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--pval_threshold"), type="numeric", default=NULL, 
                metavar="type", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--log2FC_threshold"), type="numeric", default=NULL, 
                metavar="type", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--minGS_size"), type="numeric", default=NULL, 
                metavar="type", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--maxGS_size"), type="numeric", default=NULL, 
                metavar="type", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--DE_results"), type="character", default=NULL, 
                metavar="type", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--background_genes"), type="character", default=NULL, 
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
  
  AP_function(opt)
  
  
  
}


###########################################################################

system.time( main() )