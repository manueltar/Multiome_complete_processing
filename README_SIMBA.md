# 1. To export h5ad files for SIMBA processing use:

$ bash ~/Scripts/Wraper_scripts/167_Export_RNA_and_ATAC_for_SIMBA.sh /group/soranzo/manuel.tardaguila/2025_hESC_MK_multiome/ Downstream_analysis_cluster_after_genotyping

# 2. To scan for motifs and kmers using the simba scan script (JASPAR)

$ bash ~/Scripts/Wraper_scripts/168_Simba_scan_for_kmers_motifs.sh /group/soranzo/manuel.tardaguila/2025_hESC_MK_multiome/ Downstream_analysis_cluster_after_genotyping

# 3. To run the preprocessing steps of SIMBA (QC, graph obtention, train model, compare and global embeddings) use

$ bash ~/Scripts/Wraper_scripts/170_Python_SIMBA_preprocessing.sh /group/soranzo/manuel.tardaguila/2025_hESC_MK_multiome/Downstream_analysis_cluster_after_genotyping/result_SIMBA/ QC_and_embeddings

# 4. To find GRN use the notebook

SIMBA_GRN.ipynb
