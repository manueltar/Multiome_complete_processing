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


opt = NULL

options(warn = 1)


multiVals <- function(x) paste(x,collapse=";")

heatmap_function_TFs = function(option_list)
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
  
  #### READ selected_annotations ----
  
  selected_annotations = unlist(strsplit(opt$selected_annotations, split=","))
  
  cat("selected_annotations_0\n")
  cat(sprintf(as.character(selected_annotations)))
  cat("\n")
  
  #### READ and transform annotation_file ----
  
  
  annotation_file<-read.table(file=opt$annotation_file, sep="\t", header=T)
  
  cat("annotation_file_0\n")
  cat(str(annotation_file))
  cat("\n")
  
  annotation_file_long<-unique(as.data.frame(cSplit(annotation_file,sep = '|', direction = "long",
                                           splitCols = "TF_targets"),stringsAsFactors=F))
  
  cat("annotation_file_long_0\n")
  cat(str(annotation_file_long))
  cat("\n")
  
  #### READ and transform DE_result ----
  
  
  DE_result<-readRDS(file=opt$DE_result)
  
  cat("DE_result_0\n")
  cat(str(DE_result))
  cat("\n")
  
  #### READ and transform normalised_counts ----
  
  
  normalised_counts<-readRDS(file=opt$normalised_counts)
  
  cat("normalised_counts_0\n")
  cat(str(normalised_counts))
  cat("\n")
  
  ########## define annotation selections -----------------
  
  matches <- grep(paste(selected_annotations,collapse="|"),annotation_file_long$TF_targets)
  
  toMatch<-tolower(selected_annotations)
  
  
  matches_lc <- grep(paste(toMatch,collapse="|"),annotation_file_long$TF_targets)
  
  total_matches<-unique(c(matches,matches_lc))
  
  
  if(length(total_matches) >0){
    
    
    annotation_file_long_sel<-annotation_file_long[total_matches,]
    
    cat("annotation_file_long_sel_0\n")
    cat(str(annotation_file_long_sel))
    cat("\n")
    
    #### LOOP of identities  -----
    
    
    array_identities<-unique(annotation_file_long_sel$identity)
    
    
    cat("array_identities_0\n")
    cat(str(array_identities))
    cat("\n")
    
    DEBUG<-0
    
    
    for(i in 1:length(array_identities)){
      
      identity_sel<-array_identities[i]
      
      cat("-------------------------------------------------------------------------------------------------->\t")
      cat(sprintf(as.character(identity_sel)))
      cat("\n")
      
      path_identity_sel<-paste(out,identity_sel,'/',sep='')
      
      if (file.exists(path_identity_sel)){
        
        
        
        
      }else{
        
        dir.create(path_identity_sel)
      }
      
      
      DE_result_sel<-droplevels(DE_result[which(DE_result$identity == identity_sel),])
      
      if(DEBUG == 1){
        
        cat("DE_result_sel_0\n")
        cat(str(DE_result_sel))
        cat("\n")
      }
      
      annotation_file_long_sel_identity_sel<-annotation_file_long_sel[which(annotation_file_long_sel$identity == identity_sel),]
      
      if(DEBUG == 1){
      
      cat("annotation_file_long_sel_identity_sel_0\n")
      cat(str(annotation_file_long_sel_identity_sel))
      cat("\n")
      }
      
      
      annotation_file_long_sel_identity_sel$TF_targets<-factor(annotation_file_long_sel_identity_sel$TF_targets,
                                                               levels=c(selected_annotations),
                                                               ordered = T)
      
      annotation_file_long_sel_identity_sel<-annotation_file_long_sel_identity_sel[order(annotation_file_long_sel_identity_sel$TF_targets),]
      
      if(DEBUG == 1){
        cat("annotation_file_long_sel_identity_sel_0\n")
        cat(str(annotation_file_long_sel_identity_sel))
        cat("\n")
      }
      
      annotation_file_long_sel_identity_sel.dt<-data.table(annotation_file_long_sel_identity_sel, key='gene')
      
      annotation_file_sel_collapsed<-as.data.frame(annotation_file_long_sel_identity_sel.dt[,.(TF_string=paste(TF_targets, collapse='|')), by=key(annotation_file_long_sel_identity_sel.dt)], stringsAsFactors=F)
      
      
      if(DEBUG == 1){
        cat("annotation_file_sel_collapsed_0\n")
        cat(str(annotation_file_sel_collapsed))
        cat("\n")
      }
      
      annotation_file_sel_collapsed$TF_string<-factor(annotation_file_sel_collapsed$TF_string)
      
      
      if(DEBUG == 1){
        cat("annotation_file_sel_collapsed_1\n")
        cat(str(annotation_file_sel_collapsed))
        cat("\n")
      }
      
      n_levels<-length(levels(annotation_file_sel_collapsed$TF_string))
      
      if(DEBUG == 1){
        cat("n_levels_0\n")
        cat(str(n_levels))
        cat("\n")
      }
      
      
      REP<-normalised_counts[which(normalised_counts$gene%in%annotation_file_sel_collapsed$gene &
                                     normalised_counts$identity == identity_sel),]
      
      if(DEBUG == 1){
        
        cat("REP_0\n")
        cat(str(REP))
        cat("\n")
      }
      
      
      REP_wide<-as.data.frame(pivot_wider(REP, id_cols=c('gene'), 
                                          names_from=c('identity',"clone_line"), 
                                          values_from='count',
                                          names_sep='|'), stringsAsFactors=F)
      
      if(DEBUG == 1){
        
        cat("REP_wide_0\n")
        cat(str(REP_wide))
        cat("\n")
      }
      
      GeneEXP_matrix<-as.matrix(REP_wide[,-which(colnames(REP_wide)%in%c('gene'))])
      
      row.names(GeneEXP_matrix)<-REP_wide$gene
      
      if(DEBUG == 1){
        
        cat("GeneEXP_matrix_0\n")
        cat(str(GeneEXP_matrix))
        cat("\n")
      }           
      
      annotation_col<- data.frame(matrix(vector(), length(colnames(GeneEXP_matrix)), 2,
                                         dimnames=list(c(),
                                                       c("identity","clone_line"))),stringsAsFactors=F)
      
      row.names(annotation_col)<-colnames(GeneEXP_matrix)
      
      
      if(DEBUG == 1){
        cat("annotation_col_0\n")
        cat(str(annotation_col))
        cat("\n")
      }
      
      
      annotation_col$identity<-gsub("\\|.+$","", row.names(annotation_col))
      annotation_col$clone_line<-gsub("^[^\\|]+\\|","", row.names(annotation_col))
      
      if(DEBUG == 1){
        cat("annotation_col_1\n")
        cat(str(annotation_col))
        cat("\n")
      }
      
      annotation_col$Genotype<-NA
      
      annotation_col$Genotype[which(annotation_col$clone_line%in%c('wt_1','wt_2','wt_3'))]<-'wt'
      annotation_col$Genotype[which(annotation_col$clone_line%in%c('rs62237617_1','rs62237617_2','rs62237617_3'))]<-'rs62237617'
      annotation_col$Genotype[which(annotation_col$clone_line%in%c('DNMT3A_1','DNMT3A_2','DNMT3A_3'))]<-'DNMT3A'
      annotation_col$Genotype[which(annotation_col$clone_line%in%c('rs62237617_DNMT3A_1','rs62237617_DNMT3A_2','rs62237617_DNMT3A_3'))]<-'rs62237617_DNMT3A'
      
      if(DEBUG == 1){
        cat("annotation_col_PRE\n")
        cat(str(annotation_col))
        cat("\n")
        names(summary(as.factor(annotation_col$Genotype)))
        cat("\n")
        names(summary(as.factor(annotation_col$clone_line)))
        cat("\n")
        names(summary(as.factor(annotation_col$identity)))
        cat("\n")
      }
      
      annotation_col$Genotype<-factor(annotation_col$Genotype,
                                      levels=c('wt','rs62237617','DNMT3A','rs62237617_DNMT3A'),
                                      ordered=T)
      
      annotation_col$clone_line<-factor(annotation_col$clone_line,
                                        levels=levels(normalised_counts$clone_line),
                                        ordered=T)
      
      annotation_col$identity<-factor(annotation_col$identity,
                                      levels=levels(DE_result$identity),
                                      ordered=T)
      if(DEBUG == 1){
        cat("annotation_col_POST\n")
        cat(str(annotation_col))
        cat("\n")
        names(summary(as.factor(annotation_col$Genotype)))
        cat("\n")
        names(summary(as.factor(annotation_col$clone_line)))
        cat("\n")
        names(summary(as.factor(annotation_col$identity)))
        cat("\n")
      }
      
      annotation_row = data.frame(GeneClass = annotation_file_sel_collapsed$TF_string)
      
      rownames(annotation_row) = annotation_file_sel_collapsed$gene
      
      if(DEBUG == 1){
        cat("annotation_row_0\n")
        cat(str(annotation_row))
        cat("\n")
        cat(sprintf(as.character(names(summary(as.factor(annotation_row$GeneClass))))))
        cat("\n")
        cat(sprintf(as.character(summary(as.factor(annotation_row$GeneClass)))))
        cat("\n")
      }
      
      
     
      
      
      
      vector_colors_clone_line<-c(brewer.pal(9, "Greens")[c(5)],brewer.pal(9, "Greens")[c(6)],brewer.pal(9, "Greens")[c(7)],
                                  brewer.pal(9, "Reds")[c(5)],brewer.pal(9, "Reds")[c(6)],brewer.pal(9, "Reds")[c(7)],
                                  brewer.pal(9, "Purples")[c(5)],brewer.pal(9, "Purples")[c(6)],brewer.pal(9, "Purples")[c(7)],
                                  brewer.pal(9, "Blues")[c(3)],brewer.pal(9, "Blues")[c(4)],brewer.pal(9, "Blues")[c(5)])
      
      names(vector_colors_clone_line)<-levels(annotation_col$clone_line)
      
      vector_colors_Genotype<-c(brewer.pal(9, "Greens")[c(5)],
                                brewer.pal(9, "Reds")[c(5)],
                                brewer.pal(9, "Purples")[c(5)],
                                brewer.pal(9, "Blues")[c(3)])
      
      names(vector_colors_Genotype)<-levels(annotation_col$Genotype)
      
      vector_colors_identity<-c(brewer.pal(8, "Dark2")[c(1)],
                                brewer.pal(8, "Dark2")[c(2)],
                                brewer.pal(8, "Dark2")[c(3)],
                                brewer.pal(8, "Dark2")[c(4)],
                                brewer.pal(8, "Dark2")[c(5)],
                                brewer.pal(8, "Dark2")[c(6)],
                                brewer.pal(8, "Dark2")[c(7)])
      
      
      names(vector_colors_identity)<-levels(annotation_col$identity)
      
      
      GeneClass_vector<-rep("0",n_levels)
      
      for(iteration_n_levels in 1:length(GeneClass_vector)){
        
        if(iteration_n_levels <= 8){
          
          GeneClass_vector[iteration_n_levels]<-brewer.pal(8, "Dark2")[iteration_n_levels]
          
        }else{
          
          GeneClass_vector[iteration_n_levels]<-brewer.pal(12, "Set3")[(iteration_n_levels - (8))]
          
        }#iteration_n_levels <= 8
        
        
        
      }#iteration_n_levels in 1:length(GeneClass_vector)
      
      names(GeneClass_vector)<-levels(annotation_file_sel_collapsed$TF_string)
      
      if(DEBUG == 1){
        cat("GeneClass_vector_0\n")
        cat(str(GeneClass_vector))
        cat("\n")
      }

      
      
      ann_colors <- list( clone_line = vector_colors_clone_line,
                          Genotype = vector_colors_Genotype,
                          identity =vector_colors_identity,
                          GeneClass = GeneClass_vector)
      
      
      if(DEBUG == 1){
        cat("ann_colors_0\n")
        cat(str(ann_colors))
        cat("\n")
      }
      
      
      
      
      FLAG_log_pval<-NA
      
      
      if(dim(REP_wide)[1] <= 75 & dim(REP_wide)[1] > 1){
        
        
        heatmap<-pheatmap(GeneEXP_matrix, display_numbers = FALSE,
                          show_colnames=FALSE,
                          angle_col = "0",
                          clustering_method="ward.D2",
                          fontsize_row = 8, 
                          fontsize_col = 8,
                          breaks=seq(-2,2,length.out=101),
                          color=colorRampPalette(c("blue","white","red"))(100),
                          scale="row",
                          cluster_cols=TRUE,
                          border_color='black',
                          treeheight_row=70, treeheight_col=70, cutree_cols=7,
                          annotation_col = annotation_col,
                          annotation_row = annotation_row,
                          annotation_colors = ann_colors)
        
        setwd(path_identity_sel)
        
        svgname<-paste("Heatmap_",paste(selected_annotations, collapse="_"),".svg",sep='')
        
        ggsave(svgname,plot=heatmap, device ='svg', width=13, height=13)
        
        FLAG_log_pval<-1
        
        
      }else{
        
        if(dim(REP_wide)[1] > 75){
          
          heatmap<-pheatmap(GeneEXP_matrix, display_numbers = FALSE,
                            show_colnames=FALSE,
                            show_rownames=FALSE,
                            angle_col = "0",
                            clustering_method="ward.D2",
                            fontsize_row = 8, 
                            fontsize_col = 8,
                            breaks=seq(-2,2,length.out=101),
                            color=colorRampPalette(c("blue","white","red"))(100),
                            scale="row",
                            cluster_cols=TRUE,
                            border_color='black',
                            treeheight_row=70, treeheight_col=70, cutree_cols=7,
                            annotation_col = annotation_col,
                            annotation_row = annotation_row,
                            annotation_colors = ann_colors)
          
          setwd(path_identity_sel)
          
          svgname<-paste("Heatmap_",paste(selected_annotations, collapse="_"),".svg",sep='')
          
          ggsave(svgname,plot=heatmap, device ='svg', width=13, height=13)
          
          FLAG_log_pval<-1
          
        }else{
          
          cat("no heatmap\n")
          
          FLAG_log_pval<-0
          
        }# dim(REP_wide)[1] > 50
        
      }#dim(REP_wide)[1] <= 50
      
      
      
      
      
      
      
      
      ##### logpval plot ----
      
      if(FLAG_log_pval == 1){
        selected_genes_after_heatmap_clustering<-heatmap$tree_row$labels[heatmap$tree_row$order]
        
        if(DEBUG == 1){
          cat("selected_genes_after_heatmap_clustering_0\n")
          cat(str(selected_genes_after_heatmap_clustering))
          cat("\n")
        }
        
        logpval_df<-DE_result_sel[which(DE_result_sel$gene%in%selected_genes_after_heatmap_clustering),]
        
        
        if(DEBUG == 1){
          cat("logpval_df_0\n")
          cat(str(logpval_df))
          cat("\n")
        }
        
        logpval_df$gene<-factor(logpval_df$gene, levels=rev(selected_genes_after_heatmap_clustering), ordered=T)
        
        
        
        logpval_df$SIG<-NA
        
        logpval_df$SIG[which(logpval_df$minuslog10padj >= 1.3)]<-'YES'
        logpval_df$SIG[which(logpval_df$minuslog10padj < 1.3)]<-'NO'
        
        
        logpval_df$SIG<-factor(logpval_df$SIG, levels=c('NO','YES'), ordered=T)
        
        if(DEBUG == 1){
          cat("logpval_df_1\n")
          cat(str(logpval_df))
          cat("\n")
        }
        
        
        vector_fill<-c(brewer.pal(9, "Reds")[5],brewer.pal(9, "Purples")[5],brewer.pal(9, "Blues")[3])
        
        
        logpval_dotplot<-ggplot(data=logpval_df,
                                aes(y=gene,
                                    x=contrast))+
          geom_point(aes(size=minuslog10padj, 
                         color=SIG,
                         fill=contrast),
                     stroke=1, shape=21)+
          scale_size(range = c(0,6), name='-log10pval')+
          scale_y_discrete(name=NULL)+
          scale_x_discrete(name=NULL)+
          scale_fill_manual(values=vector_fill, drop=F)+
          scale_color_manual(name='p < 0.05',values=c('gray','black'))+
          theme_classic()+
          theme(axis.title.y=element_blank(),
                axis.title.x=element_blank(),
                axis.text.y=element_text(size=8, color="black", family="sans"),
                axis.text.x=element_blank(),
                axis.line.x = element_line(size = 0.4),
                axis.ticks.x = element_blank(),
                axis.ticks.y = element_line(size = 0.4),
                axis.line.y = element_line(size = 0.4))+
          theme(legend.title = element_text(size=12),
                legend.text = element_text(size=8),
                legend.key.size = unit(0.5, 'cm'), #change legend key size
                legend.key.height = unit(0.5, 'cm'), #change legend key height
                legend.key.width = unit(0.5, 'cm'), #change legend key width
                legend.position="right")+
          ggeasy::easy_center_title()
        
        setwd(path_identity_sel)
        

        svgname<-paste("logpval_dotplot_",paste(selected_annotations, collapse="_"),".svg",sep='')
        
        ggsave(svgname,plot=logpval_dotplot, device ='svg')
        
      }#FLAG_log_pval == 1
      
      
    }# i in 1:length(array_identities)
    
    
    
    
  }#length(total_matches) >0
  
 
}

heatmap_function_other = function(option_list)
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
  
  #### READ selected_annotations_other ----
  
  selected_annotations_other = unlist(strsplit(opt$selected_annotations_other, split=","))
  
  cat("selected_annotations_other\n")
  cat(sprintf(as.character(selected_annotations_other)))
  cat("\n")
  
  #### READ and transform annotation_file ----
  
  
  annotation_file<-read.table(file=opt$annotation_file, sep="\t", header=T)
  
  cat("annotation_file_0\n")
  cat(str(annotation_file))
  cat("\n")
  
  annotation_file_long<-unique(as.data.frame(cSplit(annotation_file,sep = '|', direction = "long",
                                                    splitCols = "other"),stringsAsFactors=F))
  
  cat("annotation_file_long_0\n")
  cat(str(annotation_file_long))
  cat("\n")
  
  #### READ and transform DE_result ----
  
  
  DE_result<-readRDS(file=opt$DE_result)
  
  cat("DE_result_0\n")
  cat(str(DE_result))
  cat("\n")
  
  #### READ and transform normalised_counts ----
  
  
  normalised_counts<-readRDS(file=opt$normalised_counts)
  
  cat("normalised_counts_0\n")
  cat(str(normalised_counts))
  cat("\n")
  
  ########## define annotation selections -----------------
  
  matches <- grep(paste(selected_annotations_other,collapse="|"),annotation_file_long$other)
  
  toMatch<-tolower(selected_annotations_other)
  
  
  matches_lc <- grep(paste(toMatch,collapse="|"),annotation_file_long$other)
  
  total_matches<-unique(c(matches,matches_lc))
  
  
  if(length(total_matches) >0){
    
    
    annotation_file_long_sel<-annotation_file_long[total_matches,]
    
    cat("annotation_file_long_sel_0\n")
    cat(str(annotation_file_long_sel))
    cat("\n")
    
    #### LOOP of identities  -----
    
    
    array_identities<-unique(annotation_file_long_sel$identity)
    
    
    cat("array_identities_0\n")
    cat(str(array_identities))
    cat("\n")
    
    DEBUG<-0
    
    
    for(i in 1:length(array_identities)){
      
      identity_sel<-array_identities[i]
      
      cat("-------------------------------------------------------------------------------------------------->\t")
      cat(sprintf(as.character(identity_sel)))
      cat("\n")
      
      path_identity_sel<-paste(out,identity_sel,'/',sep='')
      
      if (file.exists(path_identity_sel)){
        
        
        
        
      }else{
        
        dir.create(path_identity_sel)
      }
      
      
      DE_result_sel<-droplevels(DE_result[which(DE_result$identity == identity_sel),])
      
      if(DEBUG == 1){
        
        cat("DE_result_sel_0\n")
        cat(str(DE_result_sel))
        cat("\n")
      }
      
      annotation_file_long_sel_identity_sel<-annotation_file_long_sel[which(annotation_file_long_sel$identity == identity_sel),]
      
      if(DEBUG == 1){
        
        cat("annotation_file_long_sel_identity_sel_0\n")
        cat(str(annotation_file_long_sel_identity_sel))
        cat("\n")
      }
      
      
      annotation_file_long_sel_identity_sel$other<-factor(annotation_file_long_sel_identity_sel$other,
                                                               levels=c(selected_annotations_other),
                                                               ordered = T)
      
      annotation_file_long_sel_identity_sel<-annotation_file_long_sel_identity_sel[order(annotation_file_long_sel_identity_sel$other),]
      
      if(DEBUG == 1){
        cat("annotation_file_long_sel_identity_sel_0\n")
        cat(str(annotation_file_long_sel_identity_sel))
        cat("\n")
      }
      
      annotation_file_long_sel_identity_sel.dt<-data.table(annotation_file_long_sel_identity_sel, key='gene')
      
      annotation_file_sel_collapsed<-as.data.frame(annotation_file_long_sel_identity_sel.dt[,.(other_string=paste(other, collapse='|')), by=key(annotation_file_long_sel_identity_sel.dt)], stringsAsFactors=F)
      
      
      if(DEBUG == 1){
        cat("annotation_file_sel_collapsed_0\n")
        cat(str(annotation_file_sel_collapsed))
        cat("\n")
      }
      
      annotation_file_sel_collapsed$other_string<-factor(annotation_file_sel_collapsed$other_string)
      
      
      if(DEBUG == 1){
        cat("annotation_file_sel_collapsed_1\n")
        cat(str(annotation_file_sel_collapsed))
        cat("\n")
      }
      
      n_levels<-length(levels(annotation_file_sel_collapsed$other_string))
      
      if(DEBUG == 1){
        cat("n_levels_0\n")
        cat(str(n_levels))
        cat("\n")
      }
      
      
      REP<-normalised_counts[which(normalised_counts$gene%in%annotation_file_sel_collapsed$gene &
                                     normalised_counts$identity == identity_sel),]
      
      if(DEBUG == 1){
        
        cat("REP_0\n")
        cat(str(REP))
        cat("\n")
      }
      
      
      REP_wide<-as.data.frame(pivot_wider(REP, id_cols=c('gene'), 
                                          names_from=c('identity',"clone_line"), 
                                          values_from='count',
                                          names_sep='|'), stringsAsFactors=F)
      
      if(DEBUG == 1){
        
        cat("REP_wide_0\n")
        cat(str(REP_wide))
        cat("\n")
      }
      
      GeneEXP_matrix<-as.matrix(REP_wide[,-which(colnames(REP_wide)%in%c('gene'))])
      
      row.names(GeneEXP_matrix)<-REP_wide$gene
      
      if(DEBUG == 1){
        
        cat("GeneEXP_matrix_0\n")
        cat(str(GeneEXP_matrix))
        cat("\n")
      }           
      
      annotation_col<- data.frame(matrix(vector(), length(colnames(GeneEXP_matrix)), 2,
                                         dimnames=list(c(),
                                                       c("identity","clone_line"))),stringsAsFactors=F)
      
      row.names(annotation_col)<-colnames(GeneEXP_matrix)
      
      
      if(DEBUG == 1){
        cat("annotation_col_0\n")
        cat(str(annotation_col))
        cat("\n")
      }
      
      
      annotation_col$identity<-gsub("\\|.+$","", row.names(annotation_col))
      annotation_col$clone_line<-gsub("^[^\\|]+\\|","", row.names(annotation_col))
      
      if(DEBUG == 1){
        cat("annotation_col_1\n")
        cat(str(annotation_col))
        cat("\n")
      }
      
      annotation_col$Genotype<-NA
      
      annotation_col$Genotype[which(annotation_col$clone_line%in%c('wt_1','wt_2','wt_3'))]<-'wt'
      annotation_col$Genotype[which(annotation_col$clone_line%in%c('rs62237617_1','rs62237617_2','rs62237617_3'))]<-'rs62237617'
      annotation_col$Genotype[which(annotation_col$clone_line%in%c('DNMT3A_1','DNMT3A_2','DNMT3A_3'))]<-'DNMT3A'
      annotation_col$Genotype[which(annotation_col$clone_line%in%c('rs62237617_DNMT3A_1','rs62237617_DNMT3A_2','rs62237617_DNMT3A_3'))]<-'rs62237617_DNMT3A'
      
      if(DEBUG == 1){
        cat("annotation_col_PRE\n")
        cat(str(annotation_col))
        cat("\n")
        names(summary(as.factor(annotation_col$Genotype)))
        cat("\n")
        names(summary(as.factor(annotation_col$clone_line)))
        cat("\n")
        names(summary(as.factor(annotation_col$identity)))
        cat("\n")
      }
      
      annotation_col$Genotype<-factor(annotation_col$Genotype,
                                      levels=c('wt','rs62237617','DNMT3A','rs62237617_DNMT3A'),
                                      ordered=T)
      
      annotation_col$clone_line<-factor(annotation_col$clone_line,
                                        levels=levels(normalised_counts$clone_line),
                                        ordered=T)
      
      annotation_col$identity<-factor(annotation_col$identity,
                                      levels=levels(DE_result$identity),
                                      ordered=T)
      if(DEBUG == 1){
        cat("annotation_col_POST\n")
        cat(str(annotation_col))
        cat("\n")
        names(summary(as.factor(annotation_col$Genotype)))
        cat("\n")
        names(summary(as.factor(annotation_col$clone_line)))
        cat("\n")
        names(summary(as.factor(annotation_col$identity)))
        cat("\n")
      }
      
      annotation_row = data.frame(GeneClass = annotation_file_sel_collapsed$other_string)
      
      rownames(annotation_row) = annotation_file_sel_collapsed$gene
      
      if(DEBUG == 1){
        cat("annotation_row_0\n")
        cat(str(annotation_row))
        cat("\n")
        cat(sprintf(as.character(names(summary(as.factor(annotation_row$GeneClass))))))
        cat("\n")
        cat(sprintf(as.character(summary(as.factor(annotation_row$GeneClass)))))
        cat("\n")
      }
      
      
      
      
      
      
      vector_colors_clone_line<-c(brewer.pal(9, "Greens")[c(5)],brewer.pal(9, "Greens")[c(6)],brewer.pal(9, "Greens")[c(7)],
                                  brewer.pal(9, "Reds")[c(5)],brewer.pal(9, "Reds")[c(6)],brewer.pal(9, "Reds")[c(7)],
                                  brewer.pal(9, "Purples")[c(5)],brewer.pal(9, "Purples")[c(6)],brewer.pal(9, "Purples")[c(7)],
                                  brewer.pal(9, "Blues")[c(3)],brewer.pal(9, "Blues")[c(4)],brewer.pal(9, "Blues")[c(5)])
      
      names(vector_colors_clone_line)<-levels(annotation_col$clone_line)
      
      vector_colors_Genotype<-c(brewer.pal(9, "Greens")[c(5)],
                                brewer.pal(9, "Reds")[c(5)],
                                brewer.pal(9, "Purples")[c(5)],
                                brewer.pal(9, "Blues")[c(3)])
      
      names(vector_colors_Genotype)<-levels(annotation_col$Genotype)
      
      vector_colors_identity<-c(brewer.pal(8, "Dark2")[c(1)],
                                brewer.pal(8, "Dark2")[c(2)],
                                brewer.pal(8, "Dark2")[c(3)],
                                brewer.pal(8, "Dark2")[c(4)],
                                brewer.pal(8, "Dark2")[c(5)],
                                brewer.pal(8, "Dark2")[c(6)],
                                brewer.pal(8, "Dark2")[c(7)])
      
      
      names(vector_colors_identity)<-levels(annotation_col$identity)
      
      
      GeneClass_vector<-rep("0",n_levels)
      
      for(iteration_n_levels in 1:length(GeneClass_vector)){
        
        if(iteration_n_levels <= 8){
          
          GeneClass_vector[iteration_n_levels]<-brewer.pal(8, "Dark2")[iteration_n_levels]
          
        }else{
          
          GeneClass_vector[iteration_n_levels]<-brewer.pal(12, "Set3")[(iteration_n_levels - (8))]
          
        }#iteration_n_levels <= 8
        
        
        
      }#iteration_n_levels in 1:length(GeneClass_vector)
      
      names(GeneClass_vector)<-levels(annotation_file_sel_collapsed$other_string)
      
      if(DEBUG == 1){
        cat("GeneClass_vector_0\n")
        cat(str(GeneClass_vector))
        cat("\n")
      }
      
      
      
      ann_colors <- list( clone_line = vector_colors_clone_line,
                          Genotype = vector_colors_Genotype,
                          identity =vector_colors_identity,
                          GeneClass = GeneClass_vector)
      
      
      if(DEBUG == 1){
        cat("ann_colors_0\n")
        cat(str(ann_colors))
        cat("\n")
      }
      
      
      
      
      FLAG_log_pval<-NA
      
      
      if(dim(REP_wide)[1] <= 75 & dim(REP_wide)[1] > 1){
        
        
        heatmap<-pheatmap(GeneEXP_matrix, display_numbers = FALSE,
                          show_colnames=FALSE,
                          angle_col = "0",
                          clustering_method="ward.D2",
                          fontsize_row = 8, 
                          fontsize_col = 8,
                          breaks=seq(-2,2,length.out=101),
                          color=colorRampPalette(c("blue","white","red"))(100),
                          scale="row",
                          cluster_cols=TRUE,
                          border_color='black',
                          treeheight_row=70, treeheight_col=70, cutree_cols=7,
                          annotation_col = annotation_col,
                          annotation_row = annotation_row,
                          annotation_colors = ann_colors)
        
        setwd(path_identity_sel)
        
        svgname<-paste("Heatmap_",paste(selected_annotations_other, collapse="_"),".svg",sep='')
        
        ggsave(svgname,plot=heatmap, device ='svg', width=13, height=13)
        
        FLAG_log_pval<-1
        
        
      }else{
        
        if(dim(REP_wide)[1] > 75){
          
          heatmap<-pheatmap(GeneEXP_matrix, display_numbers = FALSE,
                            show_colnames=FALSE,
                            show_rownames=FALSE,
                            angle_col = "0",
                            clustering_method="ward.D2",
                            fontsize_row = 8, 
                            fontsize_col = 8,
                            breaks=seq(-2,2,length.out=101),
                            color=colorRampPalette(c("blue","white","red"))(100),
                            scale="row",
                            cluster_cols=TRUE,
                            border_color='black',
                            treeheight_row=70, treeheight_col=70, cutree_cols=7,
                            annotation_col = annotation_col,
                            annotation_row = annotation_row,
                            annotation_colors = ann_colors)
          
          setwd(path_identity_sel)
          
          svgname<-paste("Heatmap_",paste(selected_annotations_other, collapse="_"),".svg",sep='')
          
          ggsave(svgname,plot=heatmap, device ='svg', width=13, height=13)
          
          FLAG_log_pval<-1
          
        }else{
          
          cat("no heatmap\n")
          
          FLAG_log_pval<-0
          
        }# dim(REP_wide)[1] > 50
        
      }#dim(REP_wide)[1] <= 50
      
      
      
      
      
      
      
      
      ##### logpval plot ----
      
      if(FLAG_log_pval == 1){
        selected_genes_after_heatmap_clustering<-heatmap$tree_row$labels[heatmap$tree_row$order]
        
        if(DEBUG == 1){
          cat("selected_genes_after_heatmap_clustering_0\n")
          cat(str(selected_genes_after_heatmap_clustering))
          cat("\n")
        }
        
        logpval_df<-DE_result_sel[which(DE_result_sel$gene%in%selected_genes_after_heatmap_clustering),]
        
        
        if(DEBUG == 1){
          cat("logpval_df_0\n")
          cat(str(logpval_df))
          cat("\n")
        }
        
        logpval_df$gene<-factor(logpval_df$gene, levels=rev(selected_genes_after_heatmap_clustering), ordered=T)
        
        
        
        logpval_df$SIG<-NA
        
        logpval_df$SIG[which(logpval_df$minuslog10padj >= 1.3)]<-'YES'
        logpval_df$SIG[which(logpval_df$minuslog10padj < 1.3)]<-'NO'
        
        
        logpval_df$SIG<-factor(logpval_df$SIG, levels=c('NO','YES'), ordered=T)
        
        if(DEBUG == 1){
          cat("logpval_df_1\n")
          cat(str(logpval_df))
          cat("\n")
        }
        
        
        vector_fill<-c(brewer.pal(9, "Reds")[5],brewer.pal(9, "Purples")[5],brewer.pal(9, "Blues")[3])
        
        
        logpval_dotplot<-ggplot(data=logpval_df,
                                aes(y=gene,
                                    x=contrast))+
          geom_point(aes(size=minuslog10padj, 
                         color=SIG,
                         fill=contrast),
                     stroke=1, shape=21)+
          scale_size(range = c(0,6), name='-log10pval')+
          scale_y_discrete(name=NULL)+
          scale_x_discrete(name=NULL)+
          scale_fill_manual(values=vector_fill, drop=F)+
          scale_color_manual(name='p < 0.05',values=c('gray','black'))+
          theme_classic()+
          theme(axis.title.y=element_blank(),
                axis.title.x=element_blank(),
                axis.text.y=element_text(size=8, color="black", family="sans"),
                axis.text.x=element_blank(),
                axis.line.x = element_line(size = 0.4),
                axis.ticks.x = element_blank(),
                axis.ticks.y = element_line(size = 0.4),
                axis.line.y = element_line(size = 0.4))+
          theme(legend.title = element_text(size=12),
                legend.text = element_text(size=8),
                legend.key.size = unit(0.5, 'cm'), #change legend key size
                legend.key.height = unit(0.5, 'cm'), #change legend key height
                legend.key.width = unit(0.5, 'cm'), #change legend key width
                legend.position="right")+
          ggeasy::easy_center_title()
        
        setwd(path_identity_sel)
        
        
        svgname<-paste("logpval_dotplot_",paste(selected_annotations_other, collapse="_"),".svg",sep='')
        
        ggsave(svgname,plot=logpval_dotplot, device ='svg')
        
      }#FLAG_log_pval == 1
      
      
    }# i in 1:length(array_identities)
    
    
    
    
  }#length(total_matches) >0
  
  
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
    make_option(c("--DE_result"), type="character", default=NULL, 
                metavar="type", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--normalised_counts"), type="character", default=NULL, 
               metavar="type", 
               help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--selected_annotations"), type="character", default=NULL, 
                metavar="type", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--selected_annotations_other"), type="character", default=NULL, 
                metavar="type", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--annotation_file"), type="character", default=NULL, 
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
  
  heatmap_function_TFs(opt)
  heatmap_function_other(opt)
  
  
  
  
}


###########################################################################

system.time( main() )