# Important. To create conda environments that are visible in jupyter notebooks:

1-Install ipkernel in environment
2-If the environment is R:
     2.1 install.packages('IRkernel')
     2.2 IRkernel::installspec(displayname = 'RGSEA')  # to register the kernel in the current R installation

# Install the conda environments in the Dependencies folder

# Recluster after genotyping and void any empty levels in the Seurat object.

$ bash ~/Scripts/Wraper_scripts/136_Recluster_after_genotyping.sh /group/soranzo/manuel.tardaguila/2025_hESC_MK_multiome/ Downstream_analysis_cluster_after_genotyping

# To explore subclustering in the new object and create the refined annotation use jupyter notebook

explore_subclustering.ipynb

# To see the DE analysis use jupyter notebook

explore_subclustering.ipynb

# To add specific transcription factor target gene sets from Dorothea run and copy them to the gmt folder (/home/manuel.tardaguila/GMT_files/msigdb_v2023.1.Hs_files_to_download_locally_ENTREZ/)

$ bash ~/Scripts/Wraper_scripts/74_ORA_exploration.sh /group/soranzo/manuel.tardaguila/Paper_bits/ ORA_exploration

# To add custom gene sets run and copy them to the gmt folder (/home/manuel.tardaguila/GMT_files/msigdb_v2023.1.Hs_files_to_download_locally_ENTREZ/)

$ bash ~/Scripts/Wraper_scripts/142_custom_genesets.sh /group/soranzo/manuel.tardaguila/gene_sets/ custom_set_for_fibrosis

# Ro run the GSEA and the ORA run:

$ bash ~/Scripts/Wraper_scripts/135_GSEA_and_ORA_for_single_cell.sh /group/soranzo/manuel.tardaguila/2025_hESC_MK_multiome/Downstream_analysis_cluster_after_genotyping/ GSEA_ORA_hESC_Day_0 /group/soranzo/manuel.tardaguila/2025_hESC_MK_multiome/Downstream_analysis_cluster_after_genotyping/DE_genes_cell_type_hESC_Day_0.tsv

# To call ATAC peaks by cell type annotation run

$ bash ~/Scripts/Wraper_scripts/138_MACS2_recall_peaks_by_cell_type.sh /group/soranzo/manuel.tardaguila/2025_hESC_MK_multiome/ Downstream_analysis_cluster_after_genotyping

# To link ATAC peaks to DE genes run

$ bash ~/Scripts/Wraper_scripts/133_Link_peaks.sh /group/soranzo/manuel.tardaguila/2025_hESC_MK_multiome/ Downstream_analysis_cluster_after_genotyping

# To do the DA anlysis see the Jupyter notebook

DA_from_new_peaks.ipynb