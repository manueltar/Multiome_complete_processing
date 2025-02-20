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

sample_array=$(echo 'MCO_01326,MCO_01327,MCO_01328,MCO_01329')


a=($(echo "$sample_array" | tr "," '\n'))

declare -a arr

for i  in "${a[@]}"
do

    sample_array_sel=$i
    echo "$sample_array_sel"


    ### Seurat_second_pass

    type=$(echo "$sample_array_sel""_""Seurat_second_pass")
    outfile_Seurat_second_pass=$(echo "$Log_files""outfile_5_""$type"".log")
    touch $outfile_Seurat_second_pass
    echo -n "" > $outfile_Seurat_second_pass
    name_Seurat_second_pass=$(echo "$type""_job")


    Rscript_Seurat_second_pass=$(echo "$Rscripts_path""405_Seurat_second_pass.R")

    sample_name=$sample_array_sel
    master_path=$MASTER_ROUTE
    preliminary_filtered=$(echo "$output_dir""$sample_array_sel""/""intermediate/""preliminary_filtered.rds")
    Uniquely_genotyped_larry_barcodes_assignments=$(echo "$MASTER_ROUTE""deconvolute_LARRY/count_and_filter/Uniquely_genotyped_larry_barcodes_assignments.rds")
    Threshold_UMIS_per_cell=$(echo '3')

    myjobid_Seurat_second_pass=$(sbatch --job-name $name_Seurat_second_pass --output=$outfile_Seurat_second_pass --partition=cpuq --time=24:00:00 --nodes=1 --ntasks-per-node=15 --mem-per-cpu=8192 --parsable --wrap="Rscript $Rscript_Seurat_second_pass --sample_name $sample_name --master_path $master_path --preliminary_filtered $preliminary_filtered --Threshold_UMIS_per_cell $Threshold_UMIS_per_cell --Uniquely_genotyped_larry_barcodes_assignments $Uniquely_genotyped_larry_barcodes_assignments --type $type --out $output_dir")
    myjobid_seff_Seurat_second_pass=$(sbatch --dependency=afterany:$myjobid_Seurat_second_pass --open-mode=append --output=$outfile_Seurat_second_pass --job-name="seff" --partition=cpuq --time=24:00:00 --nodes=1 --ntasks-per-node=1 --mem-per-cpu=128M --parsable --wrap="seff $myjobid_Seurat_second_pass >> $outfile_Seurat_second_pass")


done

conda deactivate
