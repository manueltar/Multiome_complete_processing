#!/bin/bash
 
eval "$(conda shell.bash hook)"
  

Rscripts_path=$(echo "/home/manuel.tardaguila/Scripts/R/")

MASTER_ROUTE=$1
analysis=$2

output_dir=$(echo "$MASTER_ROUTE""$analysis""/")



Log_files=$(echo "$output_dir""/""Log_files/")

conda activate multiome_NEW_downstream_analysis

##'''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''# LinkPeaks

type=$(echo "LinkPeaks")
outfile_LinkPeaks=$(echo "$Log_files""outfile_5_""$type"".log")
touch $outfile_LinkPeaks
echo -n "" > $outfile_LinkPeaks
name_LinkPeaks=$(echo "$type""_job")
seff_name=$(echo "seff""_""$type")

 

Rscript_LinkPeaks=$(echo "$Rscripts_path""424_LinkPeaks_new.R")
 
adata=$(echo "$output_dir""merged_clusters_after_genotyping_after_refined_annotation_new_peaks.rds")
DE_genes=$(echo "$output_dir""DE_genes_cell_type_hESC_Day_0.tsv")

myjobid_LinkPeaks=$(sbatch --job-name $name_LinkPeaks --output=$outfile_LinkPeaks --partition=cpuq --time=24:00:00 --nodes=1 --ntasks-per-node=8 --mem-per-cpu=8192 --parsable --wrap="Rscript $Rscript_LinkPeaks --DE_genes $DE_genes --adata $adata --type $type --out $output_dir")
myjobid_seff_LinkPeaks=$(sbatch --dependency=afterany:$myjobid_LinkPeaks --open-mode=append --output=$outfile_LinkPeaks --job-name=$seff_name --partition=cpuq --time=24:00:00 --nodes=1 --ntasks-per-node=1 --mem-per-cpu=128M --parsable --wrap="seff $myjobid_LinkPeaks >> $outfile_LinkPeaks")

conda deactivate


 
