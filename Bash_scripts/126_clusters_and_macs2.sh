#!/bin/bash

eval "$(conda shell.bash hook)"
 

Rscripts_path=$(echo "/home/manuel.tardaguila/Scripts/R/")

MASTER_ROUTE=$1
analysis=$2

output_dir=$(echo "$MASTER_ROUTE""$analysis""/")

Log_files=$(echo "$output_dir""/""Log_files/")
 
conda activate multiome_QC_DEF


### cluster_merged_object

type=$(echo "cluster_merged_object")
outfile_cluster_merged_object=$(echo "$Log_files""outfile_8_""$type"".log")
touch $outfile_cluster_merged_object
echo -n "" > $outfile_cluster_merged_object
name_cluster_merged_object=$(echo "$type""_job")


Rscript_cluster_merged_object=$(echo "$Rscripts_path""408_Clustering_of_merged_samples.R")

filtered_db_object=$(echo "$output_dir""merged_unprocessed_db_filt.rds")

myjobid_cluster_merged_object=$(sbatch --job-name $name_cluster_merged_object --output=$outfile_cluster_merged_object --partition=cpuq --time=24:00:00 --nodes=1 --ntasks-per-node=30 --mem-per-cpu=8192 --parsable --wrap="Rscript $Rscript_cluster_merged_object --filtered_db_object $filtered_db_object --type $type --out $output_dir")
myjobid_seff_cluster_merged_object=$(sbatch --dependency=afterany:$myjobid_cluster_merged_object --open-mode=append --output=$outfile_cluster_merged_object --job-name="seff" --partition=cpuq --time=24:00:00 --nodes=1 --ntasks-per-node=1 --mem-per-cpu=128M --parsable --wrap="seff $myjobid_cluster_merged_object >> $outfile_cluster_merged_object")

##########################################################################################################################

### filter_out_cluster

type=$(echo "filter_out_cluster")
outfile_filter_out_cluster=$(echo "$Log_files""outfile_9_""$type"".log")
touch $outfile_filter_out_cluster
echo -n "" > $outfile_filter_out_cluster
name_filter_out_cluster=$(echo "$type""_job")


Rscript_filter_out_cluster=$(echo "$Rscripts_path""409_Filter_out_bad_clusters.R")

db_filt_clustered=$(echo "$output_dir""merged_unprocessed_db_filt_clustered.rds")

# --dependency=afterany:

myjobid_filter_out_cluster=$(sbatch --dependency=afterany:$myjobid_cluster_merged_object --job-name $name_filter_out_cluster --output=$outfile_filter_out_cluster --partition=cpuq --time=24:00:00 --nodes=1 --ntasks-per-node=20 --mem-per-cpu=8192 --parsable --wrap="Rscript $Rscript_filter_out_cluster --db_filt_clustered $db_filt_clustered --type $type --out $output_dir")
myjobid_seff_filter_out_cluster=$(sbatch --dependency=afterany:$myjobid_filter_out_cluster --open-mode=append --output=$outfile_filter_out_cluster --job-name="seff" --partition=cpuq --time=24:00:00 --nodes=1 --ntasks-per-node=1 --mem-per-cpu=128M --parsable --wrap="seff $myjobid_filter_out_cluster >> $outfile_filter_out_cluster")

##########################################################################################################################

### MACS2_peaks_and_recluster

type=$(echo "MACS2_peaks_and_recluster")
outfile_MACS2_peaks_and_recluster=$(echo "$Log_files""outfile_10_""$type"".log")
touch $outfile_MACS2_peaks_and_recluster
echo -n "" > $outfile_MACS2_peaks_and_recluster
name_MACS2_peaks_and_recluster=$(echo "$type""_job")


Rscript_MACS2_peaks_and_recluster=$(echo "$Rscripts_path""410_Call_peaks_with_MACS2_and_reclusterize.R")

db_filt_clustered_minus_26=$(echo "$output_dir""merged_unprocessed_db_filt_clustered_minus_26.rds")
frag_file=$(echo "/group/soranzo/manuel.tardaguila/2025_hESC_MK_multiome/processing_outputs/merged.atac_fragments.tsv.gz")
db_filt_clustered_minus_26_MACS2_peaks=$(echo "$output_dir""merged_unprocessed_db_filt_clustered_minus_26_MACS2_peaks.rds")

# --dependency=afterany:$myjobid_filter_out_cluster

myjobid_MACS2_peaks_and_recluster=$(sbatch --dependency=afterany:$myjobid_filter_out_cluster --job-name $name_MACS2_peaks_and_recluster --output=$outfile_MACS2_peaks_and_recluster --partition=cpuq --time=24:00:00 --nodes=1 --ntasks-per-node=30 --mem-per-cpu=8192 --parsable --wrap="Rscript $Rscript_MACS2_peaks_and_recluster --db_filt_clustered_minus_26 $db_filt_clustered_minus_26 --frag_file $frag_file --db_filt_clustered_minus_26_MACS2_peaks $db_filt_clustered_minus_26_MACS2_peaks --type $type --out $output_dir")
myjobid_seff_MACS2_peaks_and_recluster=$(sbatch --dependency=afterany:$myjobid_MACS2_peaks_and_recluster --open-mode=append --output=$outfile_MACS2_peaks_and_recluster --job-name="seff" --partition=cpuq --time=24:00:00 --nodes=1 --ntasks-per-node=1 --mem-per-cpu=128M --parsable --wrap="seff $myjobid_MACS2_peaks_and_recluster >> $outfile_MACS2_peaks_and_recluster")

##########################################################################################################################

### Harmony_integration

type=$(echo "Harmony_integration")
outfile_Harmony_integration=$(echo "$Log_files""outfile_11_""$type"".log")
touch $outfile_Harmony_integration
echo -n "" > $outfile_Harmony_integration
name_Harmony_integration=$(echo "$type""_job")


Rscript_Harmony_integration=$(echo "$Rscripts_path""411_reclusterize_with_HARMONY.R")

db_filt_clustered_minus_26_MACS2_peaks=$(echo "$output_dir""merged_unprocessed_db_filt_clustered_minus_26_MACS2_peaks.rds")



myjobid_Harmony_integration=$(sbatch --dependency=afterany:$myjobid_MACS2_peaks_and_recluster --job-name $name_Harmony_integration --output=$outfile_Harmony_integration --partition=cpuq --time=24:00:00 --nodes=1 --ntasks-per-node=30 --mem-per-cpu=8192 --parsable --wrap="Rscript $Rscript_Harmony_integration --db_filt_clustered_minus_26_MACS2_peaks $db_filt_clustered_minus_26_MACS2_peaks --type $type --out $output_dir")
myjobid_seff_Harmony_integration=$(sbatch --dependency=afterany:$myjobid_Harmony_integration --open-mode=append --output=$outfile_Harmony_integration --job-name="seff" --partition=cpuq --time=24:00:00 --nodes=1 --ntasks-per-node=1 --mem-per-cpu=128M --parsable --wrap="seff $myjobid_Harmony_integration >> $outfile_Harmony_integration")

##########################################################################################################################

### Export_RNA_for_CellTypist

type=$(echo "Export_RNA_for_CellTypist")
outfile_Export_RNA_for_CellTypist=$(echo "$Log_files""outfile_12_""$type"".log")
touch $outfile_Export_RNA_for_CellTypist
echo -n "" > $outfile_Export_RNA_for_CellTypist
name_Export_RNA_for_CellTypist=$(echo "$type""_job")


Rscript_Export_RNA_for_CellTypist=$(echo "$Rscripts_path""413_Export_Seurat_RNA_layer_as_h5ad.R")

merged_unprocessed_db_filt_clustered_minus_26_MACS2_peaks_HARMONY_clustered=$(echo "$output_dir""merged_unprocessed_db_filt_clustered_minus_26_MACS2_peaks_HARMONY_clustered.rds")

# --dependency=afterany:$myjobid_Harmony_integration

myjobid_Export_RNA_for_CellTypist=$(sbatch --dependency=afterany:$myjobid_Harmony_integration --job-name $name_Export_RNA_for_CellTypist --output=$outfile_Export_RNA_for_CellTypist --partition=cpuq --time=24:00:00 --nodes=1 --ntasks-per-node=10 --mem-per-cpu=8192 --parsable --wrap="Rscript $Rscript_Export_RNA_for_CellTypist --merged_unprocessed_db_filt_clustered_minus_26_MACS2_peaks_HARMONY_clustered $merged_unprocessed_db_filt_clustered_minus_26_MACS2_peaks_HARMONY_clustered --type $type --out $output_dir")
myjobid_seff_Export_RNA_for_CellTypist=$(sbatch --dependency=afterany:$myjobid_Export_RNA_for_CellTypist --open-mode=append --output=$outfile_Export_RNA_for_CellTypist --job-name="seff" --partition=cpuq --time=24:00:00 --nodes=1 --ntasks-per-node=1 --mem-per-cpu=128M --parsable --wrap="seff $myjobid_Export_RNA_for_CellTypist >> $outfile_Export_RNA_for_CellTypist")


conda deactivate
 
