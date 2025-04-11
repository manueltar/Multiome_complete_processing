# To perform rpca integration see the following notebook:

RPCA_explore.ipynb

# To perform the reclustering and reach the rpca_annotation_majority_vote see notebook:

RPCA_explore.ipynb

# To run all the DE comparisons plus the GSEA and ORA run the bash script

$ bash ~/Scripts/Wraper_scripts/143_Multiome_DE_per_identity.sh /group/soranzo/manuel.tardaguila/2025_hESC_MK_multiome/Downstream_analysis_cluster_after_genotyping/AFTER_RPCA_INTEGRATION/ DE_per_identity

(For DNMT3A reference run)

$ bash ~/Scripts/Wraper_scripts/161_GSEA_and_ORA_for_single_cell_change_reference_to_DNMT3A.sh /group/soranzo/manuel.tardaguila/2025_hESC_MK_multiome/Downstream_analysis_cluster_after_genotyping/AFTER_RPCA_INTEGRATION/ DE_per_identity_change_ref_to_DNMT3A


# To run the bespoke heatmap and logpval of pathways run

$ bash ~/Scripts/Wraper_scripts/145_Multiome_bespoke_heatmaps.sh /group/soranzo/manuel.tardaguila/2025_hESC_MK_multiome/Downstream_analysis_cluster_after_genotyping/AFTER_RPCA_INTEGRATION/ DE_per_identity /group/soranzo/manuel.tardaguila/2025_hESC_MK_multiome/Downstream_analysis_cluster_after_genotyping/AFTER_RPCA_INTEGRATION/DE_per_identity/genes_GSEA_annotated.tsv

(For DNMT3A reference run)

$ bash ~/Scripts/Wraper_scripts/145_Multiome_bespoke_heatmaps.sh /group/soranzo/manuel.tardaguila/2025_hESC_MK_multiome/Downstream_analysis_cluster_after_genotyping/AFTER_RPCA_INTEGRATION/ DE_per_identity_change_ref_to_DNMT3A /group/soranzo/manuel.tardaguila/2025_hESC_MK_multiome/Downstream_analysis_cluster_after_genotyping/AFTER_RPCA_INTEGRATION/DE_per_identity_change_ref_to_DNMT3A/genes_GSEA_annotated.tsv

# To recall peaks using the the rpca_annotation run:

$ bash ~/Scripts/Wraper_scripts/138_MACS2_recall_peaks_by_cell_type_integrated_annotation.sh /group/soranzo/manuel.tardaguila/2025_hESC_MK_multiome/ Downstream_analysis_cluster_after_genotyping

# To link all DE genes to the chromatin peaks

$ bash ~/Scripts/Wraper_scripts/133_Link_peaks_NEW_peaks.sh /group/soranzo/manuel.tardaguila/2025_hESC_MK_multiome/ Downstream_analysis_cluster_after_genotyping


# To run the DA analysis:

$ bash ~/Scripts/Wraper_scripts/147_Multiome_DA_per_identity.sh /group/soranzo/manuel.tardaguila/2025_hESC_MK_multiome/Downstream_analysis_cluster_after_genotyping/AFTER_RPCA_INTEGRATION/ DA_per_identity

