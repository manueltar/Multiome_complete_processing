#!/bin/bash
 
eval "$(conda shell.bash hook)"
 
 
Rscripts_path=$(echo "/home/manuel.tardaguila/Scripts/R/")

MASTER_ROUTE=$1
analysis=$2

output_dir=$(echo "$MASTER_ROUTE""$analysis""/")

Log_files=$(echo "$output_dir""/""Log_files/")
 
conda activate multiome_QC_DEF

### recall_peaks_from_cell_type

type=$(echo "$sample_array_sel""_""recall_peaks_from_cell_type")
outfile_recall_peaks_from_cell_type=$(echo "$Log_files""outfile_4_""$type"".log")
touch $outfile_recall_peaks_from_cell_type
echo -n "" > $outfile_recall_peaks_from_cell_type
name_recall_peaks_from_cell_type=$(echo "$type""_job")


Rscript_recall_peaks_from_cell_type=$(echo "$Rscripts_path""423_MACS2_recall_peaks_from_cell_type.R")

frag_file=$(echo "/group/soranzo/manuel.tardaguila/2025_hESC_MK_multiome/processing_outputs/merged.atac_fragments.tsv.gz")
seurat_object_annotated=$(echo "/group/soranzo/manuel.tardaguila/2025_hESC_MK_multiome/Downstream_analysis_cluster_after_genotyping/merged_clusters_after_genotyping_after_refined_annotation.rds")


myjobid_recall_peaks_from_cell_type=$(sbatch --job-name $name_recall_peaks_from_cell_type --output=$outfile_recall_peaks_from_cell_type --partition=cpuq --time=24:00:00 --nodes=1 --ntasks-per-node=30 --mem-per-cpu=8192 --parsable --wrap="Rscript $Rscript_recall_peaks_from_cell_type --seurat_object_annotated $seurat_object_annotated --frag_file $frag_file --type $type --out $output_dir")
myjobid_seff_recall_peaks_from_cell_type=$(sbatch --dependency=afterany:$myjobid_recall_peaks_from_cell_type --open-mode=append --output=$outfile_recall_peaks_from_cell_type --job-name="seff" --partition=cpuq --time=24:00:00 --nodes=1 --ntasks-per-node=1 --mem-per-cpu=128M --parsable --wrap="seff $myjobid_recall_peaks_from_cell_type >> $outfile_recall_peaks_from_cell_type")


conda deactivate
