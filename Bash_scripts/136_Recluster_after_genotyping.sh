#!/bin/bash

eval "$(conda shell.bash hook)"
 
 
Rscripts_path=$(echo "/home/manuel.tardaguila/Scripts/R/")

MASTER_ROUTE=$1
analysis=$2

output_dir=$(echo "$MASTER_ROUTE""$analysis""/")

Log_files=$(echo "$output_dir""/""Log_files/")

#rm -rf $Log_files
#mkdir -p $Log_files



conda activate multiome_QC_DEF


##########################################################################################################################

### cluster_merged_object

type=$(echo "cluster_merged_object")
outfile_cluster_merged_object=$(echo "$Log_files""outfile_1_""$type"".log")
touch $outfile_cluster_merged_object
echo -n "" > $outfile_cluster_merged_object
name_cluster_merged_object=$(echo "$type""_job")


Rscript_cluster_merged_object=$(echo "$Rscripts_path""421_Clustering_after_genotyping.R")


db_genotyped=$(echo "/group/soranzo/manuel.tardaguila/2025_hESC_MK_multiome/processing_outputs/""merged_processed_presentation_updated_integrative_genotyping.rds")


myjobid_cluster_merged_object=$(sbatch --job-name $name_cluster_merged_object --output=$outfile_cluster_merged_object --partition=cpuq --time=24:00:00 --nodes=1 --ntasks-per-node=8 --mem-per-cpu=8192 --parsable --wrap="Rscript $Rscript_cluster_merged_object --db_genotyped $db_genotyped --type $type --out $output_dir")
myjobid_seff_cluster_merged_object=$(sbatch --dependency=afterany:$myjobid_cluster_merged_object --open-mode=append --output=$outfile_cluster_merged_object --job-name="seff" --partition=cpuq --time=24:00:00 --nodes=1 --ntasks-per-node=1 --mem-per-cpu=128M --parsable --wrap="seff $myjobid_cluster_merged_object >> $outfile_cluster_merged_object")

##########################################################################################################################

### harmony_genotyped

type=$(echo "harmony_genotyped")
outfile_harmony_genotyped=$(echo "$Log_files""outfile_2_""$type"".log")
touch $outfile_harmony_genotyped
echo -n "" > $outfile_harmony_genotyped
name_harmony_genotyped=$(echo "$type""_job")


Rscript_harmony_genotyped=$(echo "$Rscripts_path""422_Clust_after_geno_reclusterize_with_HARMONY.R")


merged_clusters_after_genotyping=$(echo "$output_dir""merged_clusters_after_genotyping.rds")

# --dependency=afterany:$myjobid_cluster_merged_object

myjobid_harmony_genotyped=$(sbatch --dependency=afterany:$myjobid_cluster_merged_object --job-name $name_harmony_genotyped --output=$outfile_harmony_genotyped --partition=cpuq --time=24:00:00 --nodes=1 --ntasks-per-node=8 --mem-per-cpu=8192 --parsable --wrap="Rscript $Rscript_harmony_genotyped --merged_clusters_after_genotyping $merged_clusters_after_genotyping --type $type --out $output_dir")
myjobid_seff_harmony_genotyped=$(sbatch --dependency=afterany:$myjobid_harmony_genotyped --open-mode=append --output=$outfile_harmony_genotyped --job-name="seff" --partition=cpuq --time=24:00:00 --nodes=1 --ntasks-per-node=1 --mem-per-cpu=128M --parsable --wrap="seff $myjobid_harmony_genotyped >> $outfile_harmony_genotyped")

conda deactivate


