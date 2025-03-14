# To perform rpca integration see the following notebook:

RPCA_explore.ipynb

# To perform the reclustering and reach the rpca_annotation_majority_vote see notebook:

RPCA_explore.ipynb

# To run all the DE comparisons plus the GSEA and ORA run the bash script

$ bash ~/Scripts/Wraper_scripts/143_Multiome_DE_per_identity.sh /group/soranzo/manuel.tardaguila/2025_hESC_MK_multiome/Downstream_analysis_cluster_after_genotyping/AFTER_RPCA_INTEGRATION/ DE_per_identity

# To run the bespoke heatmap and logpval of pathways run

$ bash ~/Scripts/Wraper_scripts/145_Multiome_bespoke_heatmaps.sh /group/soranzo/manuel.tardaguila/2025_hESC_MK_multiome/Downstream_analysis_cluster_after_genotyping/AFTER_RPCA_INTEGRATION/ DE_per_identity /group/soranzo/manuel.tardaguila/2025_hESC_MK_multiome/Downstream_analysis_cluster_after_genotyping/AFTER_RPCA_INTEGRATION/DE_per_identity/genes_GSEA_annotated.tsv

