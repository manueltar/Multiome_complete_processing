# General: create a conda environment from a .yml  file

$ conda env create -f renv_multiome.yml -p /home/manuel.tardaguila/conda_envs/multiome_QC_DEF

# These are the lines corresponding to the multiome mapping and Quality Control (QC)

############################################# BLOCK 1 of code ############################################
############################################# BLOCK 1 of code ############################################
############################################# BLOCK 1 of code ############################################

# 1. Map reads of the GEX and ATAC modalities using cellranger-arc count

$ sbatch ~/Scripts/sbatch/5_Cell_ranger_arccount.sh /group/soranzo/manuel.tardaguila/2025_hESC_lymph_multiome/Multiome/ <sample_id> (e.g. MCO_01330)

############################################# BLOCK 2 of code ############################################
############################################# BLOCK 2 of code ############################################
############################################# BLOCK 2 of code ############################################

## 3.1 First pass generates the initial objects (1 x run). Filters of minimum 500 RNA features, minimum 1000 ATAC features and maximum 10% Percent mithochondrial genes. Conda environment multiome_QC_DEF (see Dependencies/multiome_QC_DEF.yml).

$ bash ~/Scripts/Wraper_scripts/120_Seurat_first_v2.sh /group/soranzo/manuel.tardaguila/2025_hESC_MK_multiome/ processing_outputs

## 3.2 In parallel to 3.1 run CellBender correction of ambient RNA. IT HAS TO BE LAUNCHED FROM: processing_outputs/$SAMPLE/

$ sbatch ~/Scripts/sbatch/7_CellBender.sh /group/soranzo/manuel.tardaguila/2025_hESC_lymph_multiome/Multiome/ MCO_01330

## 3.3 In parallel to 3.1 run the python script to convert the per run ATAC peaks into 5kb matrices that can later be merged into a single object. Conda environment Manuel_ATAC (see Dependencies/Manuel_ATAC.yml).

$ bash ~/Scripts/Wraper_scripts/122_snATAC_pipeline.sh /group/soranzo/manuel.tardaguila/2025_hESC_MK_multiome/ processing_outputs

## 3.4 In parallel to 3.1 merge all the per run peak files to create a global peak reference

$ sbatch ~/Scripts/sbatch/8_merge_atac_peaks.sh /group/soranzo/manuel.tardaguila/2025_hESC_lymph_multiome/Multiome/


# 2. Map unaligned GEX reads to our reference of cell barcodes and count them. Assign barcode to a cell if 1) only one barcode association (no conflicting assignations) and 2) at least three different UMIs support the barcode assignation

  (First index the reference genome of barcodes)
  
$ bwa-mem2 index /group/soranzo/manuel.tardaguila/Multiome/RITM0023280/special_reference_files/GFP_transgene_vCHEK2_and_DNMT3A.fa

$ sbatch ~/Scripts/sbatch/6_align_to_barcodes_v2.sh /group/soranzo/manuel.tardaguila/2025_hESC_lymph_multiome/Multiome/ /group/soranzo/manuel.tardaguila/Multiome/RITM0023280/special_reference_files/GFP_transgene_vCHEK2_and_DNMT3A.fa <MCO_01330>

$ bash ~/Scripts/Wraper_scripts/119_Filter_Larry_and_graphs_v2.sh /group/soranzo/manuel.tardaguila/2025_hESC_lymph_multiome/Multiome/deconvolute_LARRY/ count_and_filter /group/soranzo/manuel.tardaguila/2025_hESC_lymph_multiome/Multiome/deconvolute_LARRY/


############################################# BLOCK 3 of code ############################################
############################################# BLOCK 3 of code ############################################
############################################# BLOCK 3 of code ############################################

## 3.5 After 3.1 has finished run Amulet to find doublets in ATAC. Conda environment multiome_QC_DEF (see Dependencies/multiome_QC_DEF.yml).

$ bash ~/Scripts/Wraper_scripts/123_Amulet.sh /group/soranzo/manuel.tardaguila/2025_hESC_MK_multiome/ processing_outputs

## 3.6 after all the previous files have been generated run the second pass per run

############################################# BLOCK 4 of code ############################################
############################################# BLOCK 4 of code ############################################
############################################# BLOCK 4 of code ############################################

$ ~/Scripts/Wraper_scripts/125_Seurat_second_pass_v2.sh

## 3.7 run the merge script to generate one object from the four separate runs

############################################# BLOCK 5 of code ############################################
############################################# BLOCK 5 of code ############################################
############################################# BLOCK 5 of code ############################################

$ bash ~/Scripts/Wraper_scripts/126_merge_pre_merged_per_sample_v2.sh /group/soranzo/manuel.tardaguila/2025_hESC_lymph_multiome/Multiome/ processing_outputs /group/soranzo/manuel.tardaguila/2025_hESC_lymph_multiome/Multiome/processing_outputs/merged.atac_fragments.tsv.gz


## 3.8 Run as a jupyter notebook: Final_QC_in_the_merged_object.ipynb . Eliminate a cluster for low quality and call final set of peaks with Macs2. Because clustering was taking a long time I run part of the jupyter notebook as a bash script

$ bash ~/Scripts/Wraper_scripts/126_clusters_and_macs2.sh /group/soranzo/manuel.tardaguila/2025_hESC_MK_multiome/ processing_outputs

############################################# targeted amplification of the GEX libraries #################################################################
############################################# targeted amplification of the GEX libraries #################################################################
############################################# targeted amplification of the GEX libraries #################################################################

# 1. Index the barcode genome with cellranger

$ cellranger mkref --fasta /group/soranzo/manuel.tardaguila/Multiome/RITM0023280/special_reference_files/GFP_transgene_vCHEK2_and_DNMT3A.fa --genes /group/soranzo/manuel.tardaguila/Multiome/RITM0023280/special_reference_files/STAR.gtf --genome GFP_transgene_vCHEK2_and_DNMT3A_cellranger

# 2. Adapt the reads to pass by input to cellranger

$ bash ~/Scripts/Wraper_scripts/127_cellranger_alignment_of_targeted_amp_GEX.sh /group/soranzo/manuel.tardaguila/2025_hESC_MK_multiome/GEX_reseq/ alignment /group/soranzo/manuel.tardaguila/Multiome/MCO_20250123/250124_A02059_0109_AHWTHYDSXC/adapter_trimmed_fastq/

# 3. Run cellranger on the adapted reads

$ sbatch ~/Scripts/Wraper_scripts/128_cellranger_MCO_1326.sh /group/soranzo/manuel.tardaguila/2025_hESC_MK_multiome/GEX_reseq/alignment/cellranger/
$ sbatch ~/Scripts/Wraper_scripts/128_cellranger_MCO_1327.sh /group/soranzo/manuel.tardaguila/2025_hESC_MK_multiome/GEX_reseq/alignment/cellranger/
$ sbatch ~/Scripts/Wraper_scripts/128_cellranger_MCO_1328.sh /group/soranzo/manuel.tardaguila/2025_hESC_MK_multiome/GEX_reseq/alignment/cellranger/
$ sbatch ~/Scripts/Wraper_scripts/128_cellranger_MCO_1329.sh /group/soranzo/manuel.tardaguila/2025_hESC_MK_multiome/GEX_reseq/alignment/cellranger/

# 4. Run CellBender to correct background of empty beads

$ sbatch ~/Scripts/Wraper_scripts/131_Cell_Bender_for_targeted_amplified_libraries.sh /group/soranzo/manuel.tardaguila/2025_hESC_MK_multiome/GEX_reseq/alignment/cellranger/ MCO_01326
$ sbatch ~/Scripts/Wraper_scripts/131_Cell_Bender_for_targeted_amplified_libraries.sh /group/soranzo/manuel.tardaguila/2025_hESC_MK_multiome/GEX_reseq/alignment/cellranger/ MCO_01327
$ sbatch ~/Scripts/Wraper_scripts/131_Cell_Bender_for_targeted_amplified_libraries.sh /group/soranzo/manuel.tardaguila/2025_hESC_MK_multiome/GEX_reseq/alignment/cellranger/ MCO_01328
$ sbatch ~/Scripts/Wraper_scripts/131_Cell_Bender_for_targeted_amplified_libraries.sh /group/soranzo/manuel.tardaguila/2025_hESC_MK_multiome/GEX_reseq/alignment/cellranger/ MCO_01329

# 5. See the jupyter notebook notebook_to_assign_barcodes.ipynb for the steps to genotype the cells
