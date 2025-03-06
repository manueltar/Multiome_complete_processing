
.libPaths()
.libPaths(new = c("/home/manuel.tardaguila/conda_envs/GSEA/lib/R/library"))
.libPaths()

suppressMessages(library(clusterProfiler))
suppressMessages(library(enrichplot))
suppressMessages(library(tidyverse))
suppressMessages(library(ggupset))
suppressMessages(library(RColorBrewer))
suppressMessages(library(pheatmap))
suppressMessages(library('org.Hs.eg.db'))
suppressMessages(library(DOSE))
suppressMessages(library("splitstackshape"))
suppressMessages(library('data.table'))
suppressMessages(library("optparse"))
suppressMessages(library('cowplot'))



opt = NULL

options(warn = 1)


multiVals <- function(x) paste(x,collapse=";")

GSEA_function = function(option_list)
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
  
  
  
  GSEA_dir<-paste0(out,'/','GSEA','/')
  
  if(file.exists(GSEA_dir)){
    
    unlink(GSEA_dir, recursive =T)
    
    dir.create(GSEA_dir)
  }else{
    
    dir.create(GSEA_dir)
  }
  
  out_path <- paste0(out,'/','GSEA','/','background_adapted/') # output path, where you want your results exported to
  
  if(file.exists(out_path)){
    
    unlink(out_path, recursive =T)
    
    dir.create(out_path)
  }else{
    
    dir.create(out_path)
  }
  
  
  cat("out_path_\n")
  cat(sprintf(as.character(out_path)))
  cat("\n")
  
  #### READ and transform path_to_GMT ----
  
  path_to_GMT = opt$path_to_GMT
  
  cat("path_to_GMT_0\n")
  cat(sprintf(as.character(path_to_GMT)))
  cat("\n")
  
  #### READ and transform search_terms ----
  
  search_terms = unlist(strsplit(opt$search_terms, split=","))
  
  cat("search_terms_\n")
  cat(sprintf(as.character(search_terms)))
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
  
  #### READ and transform DE_results ----
  
  
  DE_results<-as.data.frame(fread(file=opt$DE_results,sep="\t", header=TRUE), stringsAsFactors=F)
  
  cat("DE_results_0\n")
  cat(str(DE_results))
  cat("\n")
  cat(sprintf(as.character(names(summary(as.factor(DE_results$contrast))))))
  cat("\n")
  cat(sprintf(as.character(summary(as.factor(DE_results$contrast)))))
  cat("\n")
  
  
  ### Filter out NA padj -------------------------------
  
  DE_results_NO_NA<-DE_results[!is.na(DE_results$padj),]
  
  cat("DE_results_NO_NA_0\n")
  cat(str(DE_results_NO_NA))
  cat("\n")
  
 
  
  
  
  
  ## Get the genes that are present in your dataframe -----------------------------------------
  
  dataset_genes_df<-as.data.frame(unique(DE_results_NO_NA$gene), stringsAsFactors=F)
  colnames(dataset_genes_df)<-'gene'
  
  
  dataset_genes_df$ENTREZID <- mapIds(org.Hs.eg.db, keys=dataset_genes_df$gene, keytype="SYMBOL",
                                      column="ENTREZID", multiVals=multiVals)
  
  cat("dataset_genes_df_0\n")
  cat(str(dataset_genes_df))
  cat("\n")
  
  gmt_files <- list.files(path = path_to_GMT, pattern = '.gmt$', full.names = TRUE)
  
  cat("gmt_files_0\n")
  cat(str(gmt_files))
  cat("\n")
  
  ####### LOOP 1 save all the background files a-----------------------------------------------------------------------------------
  
  for (i in 1:length(gmt_files)){
    file <- gmt_files[i]
    
    cat("------------------------------------------->\t")
    cat(sprintf(as.character(file)))
    cat("\n")
    
    
    
    
    pwl2 <- read.gmt(file)
    
    #str(pwl2)
    
    
    
    pwl2 <- pwl2[which(pwl2$gene %in% dataset_genes_df$ENTREZID),]
    
    
    
    
    
    #str(pwl2)
    
    filename<-gsub(paste0(path_to_GMT,'/'),"",file)
    filename<-gsub("\\.gmt$","_selected.rds",filename)
    
    
    #filename <- paste('test','_',i,'.rds', sep='') #paste(gsub('c.\\.', '', gsub('.v7.5.*$', '', file)), '.RDS', sep = '')
    
    cat(sprintf(as.character(filename)))
    cat("\n")
    cat(sprintf(as.character(out_path)))
    cat("\n")
    
    
    
    
    setwd(out_path)
    
    saveRDS(pwl2, file=filename) 
    
    
    
  }#i in 1:length(gmt_files))
  
  
  # Prepare ALL genes results -----------------------------------------------
 
 
  
  DE_results_NO_NA$ENTREZID <- mapIds(org.Hs.eg.db, keys=DE_results_NO_NA$gene, keytype="SYMBOL",
                            column="ENTREZID", multiVals=multiVals)
  
  
  cat("DE_results_NO_NA_2\n")
  str(DE_results_NO_NA)
  cat("\n")
  
  
  DE_results_NO_NA_with_ENTREZ<-DE_results_NO_NA[-which(DE_results_NO_NA$ENTREZID == "NA"),]
  
  cat("DE_results_NO_NA_with_ENTREZ_0\n")
  str(DE_results_NO_NA_with_ENTREZ)
  cat("\n")
  
  
  # LOOP per contrast and cell type -----------------------------------------------------------------
  
  
  array_contrasts<-unique(DE_results_NO_NA_with_ENTREZ$contrast)
  
  cat("array_contrasts\n")
  str(array_contrasts)
  cat("\n")
  
  array_cell_types<-unique(DE_results_NO_NA_with_ENTREZ$cell_type)
  
  cat("array_cell_types\n")
  str(array_cell_types)
  cat("\n")
  
  
  DEBUG<-0
  
  bg_files <- list.files(path = out_path, pattern = '_selected.rds$', full.names = TRUE)
  
  cat("bg_files_0\n")
  str(bg_files)
  cat("\n")
  
  FINAL_results<-data.frame()
  
  List_FINAL<-list()
  
  
  for(i in 1:length(array_contrasts)){
    
    results_per_contrast<-data.frame()
    
    contrast_sel<-array_contrasts[i]
    
    cat("-------------------------------------------------------------------------------------------------------------------------------------------->\t")
    cat(sprintf(as.character(i)))
    cat("\t")
    cat(sprintf(as.character(contrast_sel)))
    cat("\n")
    
    List_cell_types<-list()
    
    for(k in 1:length(array_cell_types)){
      
      results_per_cell_type<-data.frame()
      
      cell_type_sel<-array_cell_types[k]
      
      
      cat("---------->\t")
      cat(sprintf(as.character(k)))
      cat("\t")
      cat(sprintf(as.character(cell_type_sel)))
      cat("\n")
      
      DE_results_NO_NA_with_ENTREZ_sel<-DE_results_NO_NA_with_ENTREZ[which(DE_results_NO_NA_with_ENTREZ$contrast == contrast_sel & DE_results_NO_NA_with_ENTREZ$cell_type == cell_type_sel),]
      
      if(dim(DE_results_NO_NA_with_ENTREZ_sel)[1] > 0){
        
        if(DEBUG ==1){
          
          cat("DE_results_NO_NA_with_ENTREZ_sel_0\n")
          cat(str(DE_results_NO_NA_with_ENTREZ_sel))
          cat("\n")
        }
        
        
        

        DE_results_NO_NA_with_ENTREZ_sel_ordered<-DE_results_NO_NA_with_ENTREZ_sel[order(DE_results_NO_NA_with_ENTREZ_sel$log2FoldChange, decreasing=TRUE),]
        
        if(DEBUG ==1){
          
          cat("DE_results_NO_NA_with_ENTREZ_sel_ordered_0\n")
          cat(str(DE_results_NO_NA_with_ENTREZ_sel_ordered))
          cat("\n")
        }
        
        
        List_bg_files<-list()
        
        
        for(iteration_bg_files in 1:length(bg_files)){
          
          
          
          selected_collection<-bg_files[iteration_bg_files]
          
          
          cat("--------------------------------------------------------------------->\t")
          cat(sprintf(as.character(iteration_bg_files)))
          cat("\t")
          cat(sprintf(as.character(selected_collection)))
          cat("\n")
          
          selected_collection_df<-readRDS(file=selected_collection)
          
          if(DEBUG ==1){
            
            cat("selected_collection_df_0\n")
            # str(selected_collection_df)
            cat("\n")
          }
          
          
          FLAG_custom<-length(grep(paste("Dorothea_","Custom_Soranzo_", sep="|"), selected_collection))
          
          
          
          cat("FLAG_custom_0\n")
          str(FLAG_custom)
          cat("\n")
          
          if(FLAG_custom == 0){
            
            universe<-unique(selected_collection_df$gene)
            minGSSize_spec<-10
            maxGSSize_spec<-500
            DEBUG<-0
            
            
          }else{
            
            setwd(out_path)
            
            HPO<-readRDS(file="c5.hpo.v2024.1.Hs.entrez_selected.rds") 
            
            if(DEBUG ==1){
              
              cat("HPO_0\n")
              str(HPO)
              cat("\n")
            }
            
            universe<-unique(HPO$gene)
            minGSSize_spec<-5
            maxGSSize_spec<-500
            
            DEBUG<-1
            
            if(DEBUG ==1){
              
              cat("selected_collection_df_0\n")
              str(selected_collection_df)
              cat("\n")
            }
            # universe<-unique(selected_collection_df$gene)
            
          }
          
          
          matches <- grep(paste(search_terms,collapse="|"),selected_collection_df$term)
          
          toMatch<-tolower(search_terms)
          
          
          matches_lc <- grep(paste(toMatch,collapse="|"),selected_collection_df$term)
          
          total_matches<-unique(c(matches,matches_lc))               
          
          
          if(length(total_matches) > 1){
            
            selected_collection_df_sel<-selected_collection_df[total_matches,]
            
            
            if(DEBUG ==1){
              
              cat("selected_collection_df_sel_0\n")
              str(selected_collection_df_sel)
              cat("\n")
            }
            
            #### Key function of GSEA -----------------------------------------------
            
            geneList<-DE_results_NO_NA_with_ENTREZ_sel_ordered$log2FoldChange
            
            names(geneList)<-DE_results_NO_NA_with_ENTREZ_sel_ordered$ENTREZID
            
            
            if(DEBUG ==1){
              
              cat("geneList_0\n")
              str(geneList)
              cat("\n")
            }
            
            
            
            res = GSEA(geneList, 
                       TERM2GENE=selected_collection_df_sel,                            
                       maxGSSize=maxGSSize_spec,
                       minGSSize=minGSSize_spec,
                       pvalueCutoff = 0.05,
                       pAdjustMethod = "BH")
            
            
            
            
            if(DEBUG ==1){
              
              cat("res_0\n")
              # str(res)
              cat("\n")
            }
            
            
            res_df <- res@result
            
            if(dim(res_df)[1] >0){
              
              if(DEBUG ==1){
                
                cat("res_df_0\n")
                str(res_df)
                cat("\n")
              }
              
              res_df <- res_df %>% mutate(minuslog10padj = -log10(p.adjust))
              
              if(DEBUG ==1){
                
                cat("res_df_3\n")
                str(res_df)
                cat("\n")
              }
              
              
              FLAG_ZERO<-length(which(unique(res_df$core_enrichment) == ""))
              
              
              cat("FLAG_ZERO_0\n")
              str(FLAG_ZERO)
              cat("\n")
              
              
              if(FLAG_ZERO == 0){
                
                res_df_long<-unique(as.data.frame(cSplit(res_df,sep = '/', direction = "long",
                                                         splitCols = "core_enrichment"),stringsAsFactors=F))
                
                res_df_long$core_enrichment<-as.character(res_df_long$core_enrichment)
                
                if(DEBUG ==1){
                  str(res_df_long)
                  cat("\n")
                }
                
                res_df_long$core_enrichment <- mapIds(org.Hs.eg.db, keys=res_df_long$core_enrichment, keytype="ENTREZID",
                                             column="SYMBOL", multiVals=multiVals)
                
                if(DEBUG ==1){
                  str(res_df_long)
                  cat("\n")
                }
                
                res_df_long.dt<-data.table(res_df_long, key=colnames(res_df_long)[-which(colnames(res_df_long) == 'core_enrichment')])
                
                
                res_df_long_collapsed<-as.data.frame(res_df_long.dt[,.(core_enrichment=paste(core_enrichment, collapse='/')), by=key(res_df_long.dt)], stringsAsFactors=F)
                
                if(DEBUG ==1){
                  str(res_df_long_collapsed)
                  cat("\n")
                }
                
                
                
                results_per_cell_type<-rbind(res_df_long_collapsed, results_per_cell_type)     
                
              }else{
                
                results_per_cell_type<-rbind(res_df, results_per_cell_type)          
                
                
                
              }# FLAG_ZERO == 0
              
              ##########
              
             
              
              List_bg_files[[selected_collection]]<-res
              
              names(List_bg_files)<-gsub("\\.entrez_selected\\.rds$","",gsub(paste(out_path,'/',sep=""),"",names(List_bg_files)))
              
              
              
            }#dim(res_df)[1] >0
            
            
            
            
            
          }# length(total_matches) > 1
          
          
          
          
        }#iteration_bg_files in 1:length(bg_files)
        
        
        
        if(dim(results_per_cell_type)[1] >0){
          
          results_per_cell_type$cell_type<-cell_type_sel        
          
          
          if(DEBUG ==1){
            cat("results_per_cell_type-----------------------------------------------0\n")
            str(results_per_cell_type)
            cat("\n")
            cat(sprintf(as.character(names(summary(as.factor(results_per_cell_type$cell_type))))))
            cat("\n")
            cat(sprintf(as.character(summary(as.factor(results_per_cell_type$cell_type)))))
            cat("\n")
            cat(sprintf(as.character(names(summary(as.factor(results_per_cell_type$contrast))))))
            cat("\n")
            cat(sprintf(as.character(summary(as.factor(results_per_cell_type$contrast)))))
            cat("\n")
            cat(sprintf(as.character(names(summary(as.factor(results_per_cell_type$ID))))))
            cat("\n")
            cat(sprintf(as.character(summary(as.factor(results_per_cell_type$ID)))))
            cat("\n")
          }
          
          
          results_per_contrast<-rbind(results_per_cell_type,results_per_contrast)
          
          List_cell_types[[cell_type_sel]]<-List_bg_files
          
          
        }# dim(results_per_cell_type)[1] >0    
        
        
        
      }# dim(DE_results_NO_NA_with_ENTREZ_sel)[1] > 0
      
      
      
    }# k in 1:length(array_cell_types)
    
    
    
    
    if(dim(results_per_contrast)[1] >0){             
      
      results_per_contrast$contrast<-contrast_sel
      
      if(DEBUG ==1){
        cat("results_per_contrast\n")
        str(results_per_contrast)
        cat("\n")
        cat(sprintf(as.character(names(summary(as.factor(results_per_contrast$cell_type))))))
        cat("\n")
        cat(sprintf(as.character(summary(as.factor(results_per_contrast$cell_type)))))
        cat("\n")
        cat(sprintf(as.character(names(summary(as.factor(results_per_contrast$contrast))))))
        cat("\n")
        cat(sprintf(as.character(summary(as.factor(results_per_contrast$contrast)))))
        cat("\n")
      }
      
      
      
      FINAL_results<-rbind(results_per_contrast,FINAL_results)
      
      List_FINAL[[contrast_sel]]<-List_cell_types
      
      
    }# dim(results_per_contrast)[1] >0                     
    
    
  }# i in 1:length(array_contrasts)
  
  cat("FINAL_results\n")
  str(FINAL_results)
  cat("\n")
  cat(sprintf(as.character(names(summary(as.factor(FINAL_results$cell_type))))))
  cat("\n")
  cat(sprintf(as.character(summary(as.factor(FINAL_results$cell_type)))))
  cat("\n")
  cat(sprintf(as.character(names(summary(as.factor(FINAL_results$contrast))))))
  cat("\n")
  cat(sprintf(as.character(summary(as.factor(FINAL_results$contrast)))))
  cat("\n")
  
  
  FINAL_results_SIG<-FINAL_results[which(FINAL_results$p.adjust <= pval_threshold),]
  
  FINAL_results_SIG$contrast<-factor(FINAL_results_SIG$contrast,
                                     levels=rev(array_contrasts),
                                     ordered=T)
  
  FINAL_results_SIG$cell_type<-factor(FINAL_results_SIG$cell_type,
                                     levels=rev(array_cell_types),
                                     ordered=T)
  
  
  cat("FINAL_results_SIG\n")
  str(FINAL_results_SIG)
  cat("\n")
  cat(sprintf(as.character(names(summary(as.factor(FINAL_results_SIG$cell_type))))))
  cat("\n")
  cat(sprintf(as.character(summary(as.factor(FINAL_results_SIG$cell_type)))))
  cat("\n")
  cat(sprintf(as.character(names(summary(as.factor(FINAL_results_SIG$contrast))))))
  cat("\n")
  cat(sprintf(as.character(summary(as.factor(FINAL_results_SIG$contrast)))))
  cat("\n")
  
  ##### SAVE RESULTS ----------------------------------
  
  
  setwd(GSEA_dir)
  
  write.table(FINAL_results, file="GSEA_results.tsv", sep="\t", quote=F, row.names = F)
  
  
  write.table(FINAL_results_SIG, file="GSEA_results_significant.tsv", sep="\t", quote=F, row.names = F)
  
  saveRDS(List_FINAL, file="GSEA_complete_results.rds")
  
  saveRDS(FINAL_results_SIG, file="GSEA_results_significant.rds")
  
  
  
}

Leading_edge_printer = function(option_list)
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
  
  
  
  #### READ GSEA_complete_results.rds ----

  
  GSEA_dir<-paste0(out,'/','GSEA','/')
  
  
  setwd(GSEA_dir)
  
  
  List_GSEA<-readRDS(file="GSEA_complete_results.rds")
  
  cat("List_GSEA_0\n")
  # str(List_GSEA)
  cat("\n")
  
  
  
  Leading_edge_plots_dir<-paste0(GSEA_dir,'Leading_edge_plots','/')
  
  if(file.exists(Leading_edge_plots_dir)){
    
    unlink(Leading_edge_plots_dir, recursive =T)
    
    dir.create(Leading_edge_plots_dir)
  }else{
    
    dir.create(Leading_edge_plots_dir)
  }
  
  
  cat("Leading_edge_plots_dir_0\n")
  cat(sprintf(as.character(Leading_edge_plots_dir)))
  cat("\n")
  
  
  #### iterate through the the list of lists -----
  
  
  array_contrasts<-names(List_GSEA)
  
  cat("array_contrasts_0\n")
  str(array_contrasts)
  cat("\n")
  
  
  for(i in 1:length(array_contrasts))
  {
    
    array_contrasts_sel<-array_contrasts[i]
    
    
    cat("------------------------------------->\t")
    cat(sprintf(as.character(array_contrasts_sel)))
    cat("\n")
    
    
    List_GSEA_sel<-List_GSEA[[array_contrasts_sel]]
    
    
    array_cell_types<-names(List_GSEA_sel)
    
    cat("array_cell_types_0\n")
    str(array_cell_types)
    cat("\n")
    
    for(k in 1:length(array_cell_types))
    {
      
      array_cell_types_sel<-array_cell_types[k]
      
      
      cat("----------->\t")
      cat(sprintf(as.character(array_cell_types_sel)))
      cat("\n")
      
      List_GSEA_sel_cell_type_sel<-List_GSEA_sel[[array_cell_types_sel]]
      
      array_collections<-names(List_GSEA_sel_cell_type_sel)
      
      cat("array_collections_0\n")
      str(array_collections)
      cat("\n")
      
      for(h in 1:length(array_collections)){
        
        array_collections_sel<-array_collections[h]
        
        
        cat("--->\t")
        cat(sprintf(as.character(array_collections_sel)))
        cat("\n")
        
        List_GSEA_sel_cell_type_sel_collection_sel<-List_GSEA_sel_cell_type_sel[[array_collections_sel]]
        
        array_gene_sets<-List_GSEA_sel_cell_type_sel_collection_sel$Description
        
        cat("array_gene_sets_0\n")
        str(array_gene_sets)
        cat("\n")
        
        for(l in 1:length(array_gene_sets)){
          
          array_gene_sets_sel<-array_gene_sets[l]
          
          
          cat(">\t")
          cat(sprintf(as.character(array_gene_sets_sel)))
          cat("\n")
          
          
          
          
          Leading_edge_plot<-gseaplot(List_GSEA_sel_cell_type_sel_collection_sel, geneSetID = l, title = List_GSEA_sel_cell_type_sel_collection_sel@result$Description[l])
          
          
          
          setwd(Leading_edge_plots_dir)
          
          svgname<-paste(paste(array_contrasts_sel,array_cell_types_sel,array_gene_sets_sel, sep='_'),'.svg', sep='')
          
          
          makesvg = TRUE
          
          if (makesvg == TRUE)
          {
            ggsave(svgname, plot= Leading_edge_plot,
                   device="svg")
          }
          
          
          
          
        }#l in 1:length(array_gene_sets
      }# h in 1:length(array_collections)
    }#k in 1:length(array_cell_types)
  }# i in 1:length(array_contrasts)
  
}

Lolliplot_and_gene_annotation = function(option_list)
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
  
  #### READ and transform Threshold_number_of_genes ----
  
  Threshold_number_of_genes = opt$Threshold_number_of_genes
  
  cat("Threshold_number_of_genes_\n")
  cat(sprintf(as.character(Threshold_number_of_genes)))
  cat("\n")
  
  
  
  GSEA_dir<-paste0(out,'/','GSEA','/')
  
 
  setwd(GSEA_dir)
  

  GSEA_result<-readRDS(file="GSEA_results_significant.rds")
  
  cat("GSEA_result_0\n")
  str(GSEA_result)
  cat("\n")
  
  # Maximize pvalue by ID, cell_type, contrast -----------------------------------
  
  
  GSEA_result.dt<-data.table(GSEA_result, key=c("ID","cell_type","contrast"))
  
  GSEA_result_MAX<-as.data.frame(GSEA_result.dt[,.SD[which.max(minuslog10padj)],by=key(GSEA_result.dt)], stringsAsFactors=F)
  
  cat("GSEA_result_MAX_0\n")
  str(GSEA_result_MAX)
  cat("\n")
  
  
  TF_terms<-unique(unlist(strsplit(opt$TF_terms, split=',')))
  
  cat("TF_terms_0\n")
  str(TF_terms)
  cat("\n")
  
  TF_terms_to_match<-paste(TF_terms, collapse="|")
  
  cat("TF_terms_to_match_0\n")
  str(TF_terms_to_match)
  cat("\n")
  
  
  ## Classify pathways ------------------------------------
  
  GSEA_result_MAX$PATH_class<-NA
  
  GSEA_result_MAX$PATH_class[grep(TF_terms_to_match,GSEA_result_MAX$ID)]<-'TF_targets'
  
  GSEA_result_MAX$PATH_class[-grep(TF_terms_to_match,GSEA_result_MAX$ID)]<-'other'
  
  GSEA_result_MAX$PATH_class<-factor(GSEA_result_MAX$PATH_class, levels=c('TF_targets','other'), ordered=T)
  
  
  cat("GSEA_result_MAX_1\n")
  str(GSEA_result_MAX)
  cat("\n")
  
  ## Annotate genes in pathways ------------------------------------------
  
  ind.dep<-which(colnames(GSEA_result_MAX)%in%c('core_enrichment'))
  
  
  GSEA_result_MAX.dt<-data.table(GSEA_result_MAX, key=colnames(GSEA_result_MAX)[-ind.dep])
  
  
  
  GSEA_result_MAX_with_count<-as.data.frame(GSEA_result_MAX.dt[,.(core_enrichment = core_enrichment,
                                                                  Count=length(unlist(strsplit(core_enrichment, split='/')))), by=key(GSEA_result_MAX.dt)], stringAsFactors=F)
  

  cat("GSEA_result_MAX_with_count_0\n")
  str(GSEA_result_MAX_with_count)
  cat("\n")
  
  
  
  ## Annotate genes in pathways ------------------------------------------
  
  indx<-c(which(colnames(GSEA_result_MAX_with_count)%in%c('ID','cell_type','core_enrichment','PATH_class')))
  
  GSEA_result_MAX_with_count_sub<-unique(GSEA_result_MAX_with_count[,indx])
  
  cat("GSEA_result_MAX_with_count_sub_0\n")
  str(GSEA_result_MAX_with_count_sub)
  cat("\n")
  
  GSEA_result_MAX_with_count_sub_long<-unique(as.data.frame(cSplit(GSEA_result_MAX_with_count_sub,sep = '/', direction = "long",
                                                        splitCols = "core_enrichment") , stringsAsFactors=F))
  
  colnames(GSEA_result_MAX_with_count_sub_long)[which(colnames(GSEA_result_MAX_with_count_sub_long) == 'core_enrichment')]<-'gene'
  
  cat("GSEA_result_MAX_with_count_sub_long_0\n")
  str(GSEA_result_MAX_with_count_sub_long)
  cat("\n")
  
  GSEA_result_MAX_with_count_sub_long.dt<-data.table(GSEA_result_MAX_with_count_sub_long, key=c('gene','cell_type','PATH_class'))
  
  gene_annotation<-as.data.frame(GSEA_result_MAX_with_count_sub_long.dt[,.(string_ID=paste(unique(ID), collapse='|')), by=key(GSEA_result_MAX_with_count_sub_long.dt)], stringsAsFactors=F)
  
  cat("gene_annotation_0\n")
  str(gene_annotation)
  cat("\n")
  
  gene_annotation_wide<-as.data.frame(pivot_wider(gene_annotation, id_cols=c('gene','cell_type'), names_from=PATH_class, values_from=string_ID), stringsAsFactors=F)
  
  cat("gene_annotation_wide_0\n")
  str(gene_annotation_wide)
  cat("\n")
  
  setwd(GSEA_dir)
  
  write.table(gene_annotation_wide, file="genes_GSEA_annotated.tsv", sep="\t", quote=F, row.names = F)
  
  ### Lolliplot -----------------------
  
  GSEA_result_MAX_with_count$PATH_class<-factor(as.character(GSEA_result_MAX_with_count$PATH_class), levels=rev(c('TF_targets','other')), ordered=T)
  
  
  cat("GSEA_result_MAX_with_count_REMEMBER\n")
  str(GSEA_result_MAX_with_count)
  cat("\n")
  
  
 
  
  DEBUG<-1
  
  GSEA_result_MAX_with_count_Thresholded<-GSEA_result_MAX_with_count[which(GSEA_result_MAX_with_count$Count >= Threshold_number_of_genes),]
  
  
  cat("GSEA_result_MAX_with_count_Thresholded_REMEMBER\n")
  str(GSEA_result_MAX_with_count_Thresholded)
  cat("\n")
  
  GSEA_result_MAX_with_count_Thresholded<-GSEA_result_MAX_with_count_Thresholded[order(GSEA_result_MAX_with_count_Thresholded$PATH_class),]
  
  levels_ID<-unique(as.character(GSEA_result_MAX_with_count_Thresholded$ID))
  
  GSEA_result_MAX_with_count_Thresholded$DUMMY<-factor(GSEA_result_MAX_with_count_Thresholded$ID, levels=levels_ID, ordered=T)
  
  if(DEBUG ==1){
    
    cat("GSEA_result_MAX_with_count_Thresholded_0\n")
    str(GSEA_result_MAX_with_count_Thresholded)
    cat("\n")
  }
  
  
  
  
 
  
  
  breaks_gene_sets<-as.numeric(GSEA_result_MAX_with_count_Thresholded$DUMMY)
  labels_gene_sets<-as.character(gsub("\\..+$","",GSEA_result_MAX_with_count_Thresholded$DUMMY))
  
  
  Gene_set_lolliplot<-ggplot(data=GSEA_result_MAX_with_count_Thresholded, 
                             aes(y=as.numeric(DUMMY),
                                 x=minuslog10padj)) +
    geom_segment(data=GSEA_result_MAX_with_count_Thresholded[which(GSEA_result_MAX_with_count_Thresholded$minuslog10padj > 0), ],
                 aes(y=as.numeric(DUMMY),
                     yend=as.numeric(DUMMY),
                     x=0,
                     xend=minuslog10padj,
                     color=NES),
                 size=0.8)+
    geom_point(data=GSEA_result_MAX_with_count_Thresholded[which(GSEA_result_MAX_with_count_Thresholded$minuslog10padj > 0), ],
               aes(color=NES),size=5, stroke=1, shape=21,fill="white")+
    geom_text(data=GSEA_result_MAX_with_count_Thresholded[which(GSEA_result_MAX_with_count_Thresholded$minuslog10padj > 0), ],
              aes(x=minuslog10padj, y=as.numeric(DUMMY), label=Count),color="black",size=2, family="sans",fontface="bold")+
  geom_segment(data=GSEA_result_MAX_with_count_Thresholded[which(GSEA_result_MAX_with_count_Thresholded$minuslog10padj < 0), ],
               aes(y=as.numeric(DUMMY),
                   yend=as.numeric(DUMMY),
                   x=0,
                   xend=minuslog10padj,
                   color=NES),
               size=0.8)+
    geom_point(data=GSEA_result_MAX_with_count_Thresholded[which(GSEA_result_MAX_with_count_Thresholded$minuslog10padj < 0), ],
               aes(color=NES), stroke=1, shape=21, fill="white")+
    geom_text(data=GSEA_result_MAX_with_count_Thresholded[which(GSEA_result_MAX_with_count_Thresholded$minuslog10padj < 0), ],
              aes(x=minuslog10padj, y=as.numeric(DUMMY), label=Count),color="black",size=2, family="sans",fontface="bold")+
    scale_color_gradient2(name=paste("Normalized","Enrichment","Score", sep="\n"),
                          low = "blue", high = "red",mid="white",midpoint=0,
                          na.value = NA)
  
  
  Gene_set_lolliplot <-Gene_set_lolliplot+
    theme_cowplot(font_size = 2,
                  font_family = "sans")+
    facet_grid(. ~ contrast+cell_type, scales='free_x', space='free_x', switch="y", drop=TRUE)+
    theme( strip.background = element_blank(),
           strip.placement = "outside",
           strip.text = element_text(size=5,color="black", family="sans"),
           panel.spacing = unit(0.2, "lines"),
           panel.background=element_rect(fill="white"),
           panel.border=element_rect(colour="white",size=0,5),
           panel.grid.major = element_blank(),
           panel.grid.minor = element_blank())+
    scale_x_continuous(name='-log10pval')+
    scale_y_continuous(name=NULL, breaks=breaks_gene_sets,
                       labels=labels_gene_sets)+
    theme_classic()+
    theme(axis.title=element_blank(),
          axis.title.y=element_blank(),
          axis.title.x=element_text(size=12,color="black", family="sans"),
          axis.text.y=element_text(size=10,color="black", family="sans", face='bold'),
          axis.text.x=element_text(size=10,color="black", family="sans"))+
    theme(legend.title = element_text(size=12),
          legend.text = element_text(size=10),
          legend.key.size = unit(0.5, 'cm'), #change legend key size
          legend.key.height = unit(0.5, 'cm'), #change legend key height
          legend.key.width = unit(0.5, 'cm'), #change legend key width
          legend.position="right")
  
  
  setwd(GSEA_dir)
  
  svgname<-paste(paste("Lolliplot",'GSEA',sep='_'),".svg",sep='')
  makesvg = TRUE
  
  if (makesvg == TRUE)
  {
    ggsave(svgname, plot= Gene_set_lolliplot,
           device="svg", width=13)
  }
  
  
 
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
    make_option(c("--search_terms"), type="character", default=NULL, 
                metavar="type", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--pval_threshold"), type="numeric", default=NULL, 
                metavar="type", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--log2FC_threshold"), type="numeric", default=NULL, 
                metavar="type", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--DE_results"), type="character", default=NULL, 
                metavar="type", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--TF_terms"), type="character", default=NULL, 
                metavar="type", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--Threshold_number_of_genes"), type="numeric", default=NULL, 
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
  
  GSEA_function(opt)
  Leading_edge_printer(opt)
  Lolliplot_and_gene_annotation(opt)
  
  
  
}


###########################################################################

system.time( main() )