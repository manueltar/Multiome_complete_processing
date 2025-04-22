#!/bin/bash

eval "$(conda shell.bash hook)"
 
  
Rscripts_path=$(echo "/home/manuel.tardaguila/Scripts/R/")

MASTER_ROUTE=$1
analysis=$2


##########################################################################

output_dir=$(echo "$MASTER_ROUTE""$analysis""/")

Log_files=$(echo "$output_dir""/""Log_files/")
 
conda activate multiome_QC_DEF
 
### Export_h5ad_RNA_and_ATAC

type=$(echo "Export_h5ad_RNA_and_ATAC")
outfile_Export_h5ad_RNA_and_ATAC=$(echo "$Log_files""outfile_7_""$type"".log")
touch $outfile_Export_h5ad_RNA_and_ATAC
echo -n "" > $outfile_Export_h5ad_RNA_and_ATAC
name_Export_h5ad_RNA_and_ATAC=$(echo "$type""_job")


Rscript_Export_h5ad_RNA_and_ATAC=$(echo "$Rscripts_path""461_export_h5ad_RNA_and_ATAC.R")


Seurat_object=$(echo "/group/soranzo/manuel.tardaguila/2025_hESC_MK_multiome/Downstream_analysis_cluster_after_genotyping/""merged_clusters_after_genotyping_after_refined_annotation_rpca_new_peaks.rds")

myjobid_Export_h5ad_RNA_and_ATAC=$(sbatch --job-name $name_Export_h5ad_RNA_and_ATAC --output=$outfile_Export_h5ad_RNA_and_ATAC --partition=cpuq --time=24:00:00 --nodes=1 --ntasks-per-node=20 --mem-per-cpu=4096 --parsable --wrap="Rscript $Rscript_Export_h5ad_RNA_and_ATAC --Seurat_object $Seurat_object --type $type --out $output_dir")
myjobid_seff_Export_h5ad_RNA_and_ATAC=$(sbatch --dependency=afterany:$myjobid_Export_h5ad_RNA_and_ATAC --open-mode=append --output=$outfile_Export_h5ad_RNA_and_ATAC --job-name="seff" --partition=cpuq --time=24:00:00 --nodes=1 --ntasks-per-node=1 --mem-per-cpu=128M --parsable --wrap="seff $myjobid_Export_h5ad_RNA_and_ATAC >> $outfile_Export_h5ad_RNA_and_ATAC")

conda deactivate
