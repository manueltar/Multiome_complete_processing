#!/bin/bash
 
eval "$(conda shell.bash hook)"
 
 
Rscripts_path=$(echo "/home/manuel.tardaguila/Scripts/R/")

MASTER_ROUTE=$1
analysis=$2

output_dir=$(echo "$MASTER_ROUTE""$analysis""/")

Log_files=$(echo "$output_dir""/""Log_files/")
 
conda activate multiome_QC_DEF

### RPCA_integration

type=$(echo "RPCA_integration")
outfile_RPCA_integration=$(echo "$Log_files""outfile_9_""$type"".log")
touch $outfile_RPCA_integration
echo -n "" > $outfile_RPCA_integration
name_RPCA_integration=$(echo "$type""_job")


Rscript_RPCA_integration=$(echo "$Rscripts_path""456_rpca_integration.R")

filt_clustered_QCed_cell_annotated=$(echo "/group/soranzo/manuel.tardaguila/2025_hESC_lymph_multiome/Multiome/processing_outputs/merged_unprocessed_db_filt_clustered_QCed_cell_annotated.rds")

# --ntasks-per-node=30 --mem-per-cpu=8192

myjobid_RPCA_integration=$(sbatch --job-name $name_RPCA_integration --output=$outfile_RPCA_integration --partition=cpuq --time=24:00:00 --nodes=1 --ntasks-per-node=30 --mem-per-cpu=12000 --parsable --wrap="Rscript $Rscript_RPCA_integration --filt_clustered_QCed_cell_annotated $filt_clustered_QCed_cell_annotated --type $type --out $output_dir")
myjobid_seff_RPCA_integration=$(sbatch --dependency=afterany:$myjobid_RPCA_integration --open-mode=append --output=$outfile_RPCA_integration --job-name="seff" --partition=cpuq --time=24:00:00 --nodes=1 --ntasks-per-node=1 --mem-per-cpu=128M --parsable --wrap="seff $myjobid_RPCA_integration >> $outfile_RPCA_integration")

### new_clusters_and_marker_genes

type=$(echo "new_clusters_and_marker_genes")
outfile_new_clusters_and_marker_genes=$(echo "$Log_files""outfile_10_""$type"".log")
touch $outfile_new_clusters_and_marker_genes
echo -n "" > $outfile_new_clusters_and_marker_genes
name_new_clusters_and_marker_genes=$(echo "$type""_job")


Rscript_new_clusters_and_marker_genes=$(echo "$Rscripts_path""459_rpca_integration_partII_new_clusters_and_marker_genes.R")

db_filt_clustered_QCed_cell_annotated_rpca_integrate=$(echo "$output_dir""merged_unprocessed_db_filt_clustered_QCed_cell_annotated_rpca_integrate.rds")

#  --dependency=afterany:$myjobid_RPCA_integration

myjobid_new_clusters_and_marker_genes=$(sbatch --dependency=afterany:$myjobid_RPCA_integration --job-name $name_new_clusters_and_marker_genes --output=$outfile_new_clusters_and_marker_genes --partition=cpuq --time=72:00:00 --nodes=1 --ntasks-per-node=20 --mem-per-cpu=8192 --parsable --wrap="Rscript $Rscript_new_clusters_and_marker_genes --db_filt_clustered_QCed_cell_annotated_rpca_integrate $db_filt_clustered_QCed_cell_annotated_rpca_integrate --type $type --out $output_dir")
myjobid_seff_new_clusters_and_marker_genes=$(sbatch --dependency=afterany:$myjobid_new_clusters_and_marker_genes --open-mode=append --output=$outfile_new_clusters_and_marker_genes --job-name="seff" --partition=cpuq --time=24:00:00 --nodes=1 --ntasks-per-node=1 --mem-per-cpu=128M --parsable --wrap="seff $myjobid_new_clusters_and_marker_genes >> $outfile_new_clusters_and_marker_genes")

 
conda deactivate
