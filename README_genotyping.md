############################################# genotyping directly from GEX libraries #################################################################
############################################# genotyping directly from GEX libraries #################################################################
############################################# genotyping directly from GEX libraries #################################################################

0. Index the reference genome of barcodes

$ bwa-mem2 index /group/soranzo/manuel.tardaguila/Multiome/RITM0023280/special_reference_files/GFP_transgene_vCHEK2_and_DNMT3A.fa

1. First align unaligned reads with bwa-mem

	sbatch ~/Scripts/sbatch/6_align_to_barcodes_v2.sh /group/soranzo/manuel.tardaguila/2025_hESC_lymph_multiome/Multiome/ /group/soranzo/manuel.tardaguila/Multiome/RITM0023280/special_reference_files/GFP_transgene_vCHEK2_and_DNMT3A.fa <sample_id>

2. Then process the aligned barcodes to get the asignation at different thresholds ofr minimum of barcodes mapped

$ bash ~/Scripts/Wraper_scripts/119_Filter_Larry_and_graphs.sh /group/soranzo/manuel.tardaguila/2025_hESC_MK_multiome/deconvolute_LARRY/ count_and_filter

or

$ bash ~/Scripts/Wraper_scripts/119_Filter_Larry_and_graphs_v2.sh /group/soranzo/manuel.tardaguila/2025_hESC_lymph_multiome/Multiome/deconvolute_LARRY/ count_and_filter /group/soranzo/manuel.tardaguila/2025_hESC_lymph_multiome/Multiome/deconvolute_LARRY/


############################################# targeted amplification of the GEX libraries #################################################################
############################################# targeted amplification of the GEX libraries #################################################################
############################################# targeted amplification of the GEX libraries #################################################################

# 1. Index the barcode genome with cellranger

$ cellranger mkref --fasta /group/soranzo/manuel.tardaguila/Multiome/RITM0023280/special_reference_files/GFP_transgene_vCHEK2_and_DNMT3A.fa --genes /group/soranzo/manuel.tardaguila/Multiome/RITM0023280/special_reference_files/STAR.gtf --genome GFP_transgene_vCHEK2_and_DNMT3A_cellranger

# 2. Adapt the reads to pass by input to cellranger

$ bash ~/Scripts/Wraper_scripts/127_cellranger_alignment_of_targeted_amp_GEX.sh /group/soranzo/manuel.tardaguila/2025_hESC_MK_multiome/GEX_reseq/ alignment /group/soranzo/manuel.tardaguila/Multiome/MCO_20250123/250124_A02059_0109_AHWTHYDSXC/adapter_trimmed_fastq/

$ bash ~/Scripts/Wraper_scripts/146_cellranger_alignment_of_targeted_amp_GEX_lymph.sh /group/soranzo/manuel.tardaguila/2025_hESC_lymph_multiome/ targeted_am\
plicon_GEX /group/soranzo/manuel.tardaguila/Multiome/QC_20250310/250310_A02059_0119_AH2N2TDMX2/adapter_trimmed_fastq/

$ bash ~/Scripts/Wraper_scripts/171_cellranger_alignment_of_targeted_amp_GEX_lymph.sh /group/soranzo/manuel.tardaguila/2025_hESC_lymph_multiome/ targeted_am\
plicon_GEX /group/soranzo/manuel.tardaguila/Multiome/RITM0029357/250414_A01481_0283_BHMKVYDRX5/fastq_raw/



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
