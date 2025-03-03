
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

ORA_function = function(option_list)
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
  
  
  
  ORA_dir<-paste0(out,'/','ORA','/')
  
  if(file.exists(ORA_dir)){
    
    unlink(ORA_dir, recursive =T)
    
    dir.create(ORA_dir)
  }else{
    
    dir.create(ORA_dir)
  }
  
  out_path <- paste0(out,'/','ORA','/','background_adapted/') # output path, where you want your results exported to
  
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
  
  ############ create  diffexpressed column -----------------
  
  
  DE_results_NO_NA <- DE_results_NO_NA %>% mutate(diffexpressed = case_when(
    log2FoldChange > log2FC_threshold & padj < pval_threshold ~ 'UP',
    log2FoldChange < log2FC_threshold & padj < pval_threshold ~ 'DOWN',
    padj > 0.05 ~ 'NO'
  ))
  
  cat("DE_results_NO_NA_1\n")
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
  
  
  # Prepare deg results -----------------------------------------------
 
  DE_results_NO_NA_SIG<-DE_results_NO_NA[which(DE_results_NO_NA$diffexpressed !='NO'), ]
  
  cat("DE_results_NO_NA_SIG_1\n")
  str(DE_results_NO_NA_SIG)
  cat("\n")
  
  DE_results_NO_NA_SIG$ENTREZID <- mapIds(org.Hs.eg.db, keys=DE_results_NO_NA_SIG$gene, keytype="SYMBOL",
                            column="ENTREZID", multiVals=multiVals)
  
  cat("DE_results_NO_NA_SIG_2\n")
  str(DE_results_NO_NA_SIG)
  cat("\n")
  
  
  # LOOP per contrast and cell type -----------------------------------------------------------------
  
  
  array_contrasts<-unique(DE_results_NO_NA_SIG$contrast)
  
  cat("array_contrasts\n")
  str(array_contrasts)
  cat("\n")
  
  array_cell_types<-unique(DE_results_NO_NA_SIG$cell_type)
  
  cat("array_cell_types\n")
  str(array_cell_types)
  cat("\n")
  
  
  DEBUG<-0
  
  bg_files <- list.files(path = out_path, pattern = '_selected.rds$', full.names = TRUE)
  
  cat("bg_files_0\n")
  str(bg_files)
  cat("\n")
  
  FINAL_results<-data.frame()
  
  
  #### key option in the ORA https://github.com/YuLab-SMU/clusterProfiler/issues/283 -------------------------------------
  
  options(enrichment_force_universe=TRUE)
  
  
  for(i in 1:length(array_contrasts)){
    
    results_per_contrast<-data.frame()
    
    contrast_sel<-array_contrasts[i]
    
    cat("-------------------------------------------------------------------------------------------------------------------------------------------->\t")
    cat(sprintf(as.character(i)))
    cat("\t")
    cat(sprintf(as.character(contrast_sel)))
    cat("\n")
    
    for(k in 1:length(array_cell_types)){
      
      results_per_cell_type<-data.frame()
      
      cell_type_sel<-array_cell_types[k]
      
      
      cat("---------->\t")
      cat(sprintf(as.character(k)))
      cat("\t")
      cat(sprintf(as.character(cell_type_sel)))
      cat("\n")
      
      DE_results_NO_NA_SIG_sel<-DE_results_NO_NA_SIG[which(DE_results_NO_NA_SIG$contrast == contrast_sel & DE_results_NO_NA_SIG$cell_type == cell_type_sel),]
      
      if(dim(DE_results_NO_NA_SIG_sel)[1] > 0){
        
        if(DEBUG ==1){
          
          cat("DE_results_NO_NA_SIG_sel_0\n")
          cat(str(DE_results_NO_NA_SIG_sel))
          cat("\n")
        }
        
        # Split the dataframe into a list of sub-dataframes: upregulated, downregulated genes
        
        deg_results_list <- split(DE_results_NO_NA_SIG_sel, DE_results_NO_NA_SIG_sel$diffexpressed)
        
        if(DEBUG ==1){
          
          cat("deg_results_list_0\n")
          cat(str(deg_results_list))
          cat("\n")
        }
        
        
        
        
        
        for(iteration_bg_files in 1:length(bg_files)){
          
          
          
          selected_collection<-bg_files[iteration_bg_files]
          
          
          cat("--->\t")
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
          
          
          FLAG_Dorothea<-length((grep("Dorothea_", selected_collection)))
          
          if(DEBUG ==1){
            
            cat("FLAG_Dorothea_0\n")
            str(FLAG_Dorothea)
            cat("\n")
          }
          
          if(FLAG_Dorothea == 0){
            
            universe<-unique(selected_collection_df$gene)
            minGSSize_spec<-10
            maxGSSize_spec<-500
            
          }else{
            
            setwd(out_path)
            
            HPO<-readRDS(file="c5.hpo.v2024.1.Hs.entrez_selected.rds") 
            
            if(DEBUG ==1){
              
              cat("HPO_0\n")
              str(HPO)
              cat("\n")
            }
            
            universe<-unique(HPO$gene)
            minGSSize_spec<-10
            maxGSSize_spec<-500
            # universe<-unique(selected_collection_df$gene)
            
          }
          
          
          matches <- grep(paste(search_terms,collapse="|"),selected_collection_df$term)
          
          toMatch<-tolower(search_terms)
          
          
          matches_lc <- grep(paste(toMatch,collapse="|"),selected_collection_df$term)
          
          total_matches<-unique(c(matches,matches_lc))               
          
          
          if(length(total_matches) > 1){
            
            selected_collection_df_sel<-selected_collection_df[total_matches,]
            
            #### Key function of ORA -----------------------------------------------
            
            
            res<-lapply(names(deg_results_list),
                        function(x) enricher(gene= deg_results_list[[x]]$ENTREZID, 
                                             TERM2GENE = selected_collection_df_sel, 
                                             universe=universe,
                                             maxGSSize=maxGSSize_spec,
                                             minGSSize=minGSSize_spec,
                                             pvalueCutoff = 0.05,
                                             pAdjustMethod = "BH",
                                             qvalueCutoff = 0.25))
            
            
            
            names(res) <- names(deg_results_list)
            
            if(DEBUG ==1){
              
              cat("res_0\n")
              #str(res)
              cat("\n")
            }
            
            ### Filter null results
            
            filtered<-Filter(Negate(is.null), res)
            
            if(DEBUG ==1){
              cat("filtered_0\n")
              #str(filtered)
              cat("\n")
            }
            
            if(length(filtered) > 0){
              
              res_df <- lapply(names(filtered), 
                               function(x) rbind(filtered[[x]]@result))
              
              if(DEBUG ==1){
                
                cat("res_df_0\n")
                str(res_df)
                cat("\n")
              }
              
              names(res_df) <- names(filtered)
              
              if(DEBUG ==1){
                
                cat("res_df_1\n")
                str(res_df)
                cat("\n")
              }
              
              res_df <- do.call(rbind, res_df)
              
              if(DEBUG ==1){
                
                cat("res_df_2\n")
                str(res_df)
                cat("\n")
              }
              
              res_df <- res_df %>% mutate(minuslog10padj = -log10(p.adjust),
                                          diffexpressed = gsub('\\..+$', '', rownames(res_df)))                    
              
              if(DEBUG ==1){
                
                cat("res_df_3\n")
                str(res_df)
                cat("\n")
              }
              
              
              
              res_df_long<-unique(as.data.frame(cSplit(res_df,sep = '/', direction = "long",
                                                       splitCols = "geneID"),stringsAsFactors=F))
              
              if(dim(res_df_long)[1] >0){
                
                res_df_long$geneID<-as.character(res_df_long$geneID)
                
                if(DEBUG ==1){
                  str(res_df_long)
                  cat("\n")
                }
                
                res_df_long$geneID <- mapIds(org.Hs.eg.db, keys=res_df_long$geneID, keytype="ENTREZID",
                                             column="SYMBOL", multiVals=multiVals)
                
                if(DEBUG ==1){
                  str(res_df_long)
                  cat("\n")
                }
                
                res_df_long.dt<-data.table(res_df_long, key=colnames(res_df_long)[-which(colnames(res_df_long) == 'geneID')])
                
                
                res_df_long_collapsed<-as.data.frame(res_df_long.dt[,.(geneID=paste(geneID, collapse='/')), by=key(res_df_long.dt)], stringsAsFactors=F)
                
                if(DEBUG ==1){
                  str(res_df_long_collapsed)
                  cat("\n")
                }
                
                
                
                results_per_cell_type<-rbind(res_df_long_collapsed, results_per_cell_type)          
                
              }else{
                
                results_per_cell_type<-rbind(res_df, results_per_cell_type)          
                
              }# dim(res_df_long)[1] >0
              
                 
              
              
            }# length(filtered) > 0     
            
            
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
            cat(sprintf(as.character(names(summary(as.factor(results_per_cell_type$diffexpressed))))))
            cat("\n")
            cat(sprintf(as.character(summary(as.factor(results_per_cell_type$diffexpressed)))))
            cat("\n")
            cat(sprintf(as.character(names(summary(as.factor(results_per_cell_type$ID))))))
            cat("\n")
            cat(sprintf(as.character(summary(as.factor(results_per_cell_type$ID)))))
            cat("\n")
          }
          
          
          results_per_contrast<-rbind(results_per_cell_type,results_per_contrast)
          
        }# dim(results_per_cell_type)[1] >0    
        
        
        
      }# dim(DE_results_NO_NA_SIG_sel)[1] > 0
      
      
      
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
  
  
  setwd(ORA_dir)
  
  write.table(FINAL_results, file="ORA_results.tsv", sep="\t", quote=F, row.names = F)
  
  
  write.table(FINAL_results_SIG, file="ORA_results_significant.tsv", sep="\t", quote=F, row.names = F)
  
 
  
  
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
  
  
  
  ORA_dir<-paste0(out,'/','ORA','/')
  
 
  setwd(ORA_dir)
  

  ORA_result<-read.table(file="ORA_results_significant.tsv", sep="\t", header=TRUE)
  
  cat("ORA_result_0\n")
  str(ORA_result)
  cat("\n")
  
  # Maximize pvalue by ID, diffexpressed, cell_type, contrast -----------------------------------
  
  
  ORA_result.dt<-data.table(ORA_result, key=c("ID","diffexpressed","cell_type","contrast"))
  
  ORA_result_MAX<-as.data.frame(ORA_result.dt[,.SD[which.max(minuslog10padj)],by=key(ORA_result.dt)], stringsAsFactors=F)
  
  cat("ORA_result_MAX_0\n")
  str(ORA_result_MAX)
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
  
  ORA_result_MAX$PATH_class<-NA
  
  ORA_result_MAX$PATH_class[grep(TF_terms_to_match,ORA_result_MAX$ID)]<-'TF_targets'
  
  ORA_result_MAX$PATH_class[-grep(TF_terms_to_match,ORA_result_MAX$ID)]<-'other'
  
  ORA_result_MAX$PATH_class<-factor(ORA_result_MAX$PATH_class, levels=c('TF_targets','other'), ordered=T)
  
  
  cat("ORA_result_MAX_1\n")
  str(ORA_result_MAX)
  cat("\n")
  
  ## Annotate genes in pathways ------------------------------------------
  
  indx<-c(which(colnames(ORA_result_MAX)%in%c('ID','diffexpressed','cell_type','geneID','PATH_class')))
  
  ORA_result_MAX_sub<-unique(ORA_result_MAX[,indx])
  
  cat("ORA_result_MAX_sub_0\n")
  str(ORA_result_MAX_sub)
  cat("\n")
  
  ORA_result_MAX_sub_long<-unique(as.data.frame(cSplit(ORA_result_MAX_sub,sep = '/', direction = "long",
                                                        splitCols = "geneID") , stringsAsFactors=F))
  
  colnames(ORA_result_MAX_sub_long)[which(colnames(ORA_result_MAX_sub_long) == 'geneID')]<-'gene'
  
  cat("ORA_result_MAX_sub_long_0\n")
  str(ORA_result_MAX_sub_long)
  cat("\n")
  
  ORA_result_MAX_sub_long.dt<-data.table(ORA_result_MAX_sub_long, key=c('gene','cell_type','diffexpressed','PATH_class'))
  
  gene_annotation<-as.data.frame(ORA_result_MAX_sub_long.dt[,.(string_ID=paste(unique(ID), collapse='|')), by=key(ORA_result_MAX_sub_long.dt)], stringsAsFactors=F)
  
  cat("gene_annotation_0\n")
  str(gene_annotation)
  cat("\n")
  
  gene_annotation_wide<-as.data.frame(pivot_wider(gene_annotation, id_cols=c('gene','cell_type','diffexpressed'), names_from=PATH_class, values_from=string_ID), stringsAsFactors=F)
  
  cat("gene_annotation_wide_0\n")
  str(gene_annotation_wide)
  cat("\n")
  
  setwd(ORA_dir)
  
  write.table(gene_annotation_wide, file="genes_ORA_annotated.tsv", sep="\t", quote=F, row.names = F)
  
  ### Lolliplot -----------------------
  
  ORA_result_MAX$PATH_class<-factor(as.character(ORA_result_MAX$PATH_class), levels=rev(c('TF_targets','other')), ordered=T)
  
  
  cat("ORA_result_MAX_REMEMBER\n")
  str(ORA_result_MAX)
  cat("\n")
  
  
  array_diffexpressed<-unique(ORA_result_MAX$diffexpressed)
  
  cat("array_diffexpressed_0\n")
  str(array_diffexpressed)
  cat("\n")
  
  
  
  DEBUG<-1
  
  
  for(i in 1:length(array_diffexpressed)){
    
    array_diffexpressed_sel<-array_diffexpressed[i]
    
    cat("--------------------------->\t")
    cat(sprintf(array_diffexpressed_sel))
    cat("\n")
    
    ORA_result_MAX_sel<-ORA_result_MAX[which(ORA_result_MAX$diffexpressed == array_diffexpressed_sel & ORA_result_MAX$Count >= Threshold_number_of_genes),]
    
    ORA_result_MAX_sel<-ORA_result_MAX_sel[order(ORA_result_MAX_sel$PATH_class),]
    
    levels_ID<-unique(as.character(ORA_result_MAX_sel$ID))
    
    ORA_result_MAX_sel$DUMMY<-factor(ORA_result_MAX_sel$ID, levels=levels_ID, ordered=T)
    
    if(DEBUG ==1){
      
      cat("ORA_result_MAX_sel_0\n")
      str(ORA_result_MAX_sel)
      cat("\n")
    }
    
    color_selected<-NA
    
    if(array_diffexpressed_sel == 'UP'){
      
      color_selected<-'red'
      
    }else{
      
      if(array_diffexpressed_sel == 'DOWN'){
        
        color_selected<-'blue'
        
      }
    }#array_diffexpressed_sel == 'UP'
    
    if(DEBUG ==1){
      
      cat("color_selected_0\n")
      str(color_selected)
      cat("\n")
    }
    
    breaks_gene_sets<-as.numeric(ORA_result_MAX_sel$DUMMY)
    labels_gene_sets<-as.character(gsub("\\..+$","",ORA_result_MAX_sel$DUMMY))
    
    
    Gene_set_lolliplot<-ggplot(data=ORA_result_MAX_sel, 
                               aes(y=as.numeric(DUMMY),
                                   x=minuslog10padj)) +
      geom_segment(data=ORA_result_MAX_sel,
                   aes(y=as.numeric(DUMMY),
                       yend=as.numeric(DUMMY),
                       x=0,
                       xend=minuslog10padj),
                   color=color_selected,
                   size=0.8)+
      geom_point(size=5, stroke=1, shape=21, color=color_selected, fill="white")+
      geom_text(data=ORA_result_MAX_sel,
                aes(x=minuslog10padj, y=as.numeric(DUMMY), label=Count),color="black",size=2, family="sans",fontface="bold")
    
    
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
            axis.title.x=element_text(size=8,color="black", family="sans"),
            axis.text.y=element_text(size=6,color="black", family="sans", face='bold'),
            axis.text.x=element_text(size=6,color="black", family="sans"),
            axis.line.x = element_line(size = 0.4),
            axis.ticks.x = element_line(size = 0.4),
            axis.ticks.y = element_line(size = 0.4),
            axis.line.y = element_line(size = 0.4))+
      theme(legend.title = element_blank(),
            legend.text = element_text(size=6),
            legend.key.size = unit(0.5, 'cm'), #change legend key size
            legend.key.height = unit(0.5, 'cm'), #change legend key height
            legend.key.width = unit(0.5, 'cm'), #change legend key width
            legend.position="bottom")+
      guides(fill=guide_legend(nrow=1,byrow=TRUE))
    
    
    setwd(ORA_dir)
    
    svgname<-paste(paste("Lolliplot",array_diffexpressed_sel,sep='_'),".svg",sep='')
    makesvg = TRUE
    
    if (makesvg == TRUE)
    {
      ggsave(svgname, plot= Gene_set_lolliplot,
             device="svg")
    }
    
    
    
  }#for(i in 1:length(array_diffexpressed))
  
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
  
  ORA_function(opt)
  Lolliplot_and_gene_annotation(opt)
  
  
  
}


###########################################################################

system.time( main() )