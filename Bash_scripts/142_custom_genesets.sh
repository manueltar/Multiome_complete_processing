#!/bin/bash
 
MASTER_ROUTE=$1
analysis=$2


Rscripts_path=$(echo "/home/manuel.tardaguila/Scripts/R/")
module load R/4.1.0


bashrc_file=$(echo "/home/manuel.tardaguila/.bashrc")

source $bashrc_file
eval "$(conda shell.bash hook)"


output_dir=$(echo "$MASTER_ROUTE""$analysis""/")

Log_files=$(echo "$output_dir""/""Log_files/")

rm -rf $Log_files
mkdir -p $Log_files

type=$(echo "Custom_print_gmt""_""$analysis")
outfile_Custom_print_gmt=$(echo "$Log_files""outfile_1_""$type"".log")
touch $outfile_Custom_print_gmt
echo -n "" > $outfile_Custom_print_gmt
name_Custom_print_gmt=$(echo "$type""_job")
seff_name=$(echo "seff""_""$type")

Rscript_Custom_print_gmt=$(echo "$Rscripts_path""432_custom_gene_sets_create_gmt_v2.R")


Table_of_gene_sets=$(echo "/group/soranzo/manuel.tardaguila/gene_sets/custom_gene_sets.txt")


myjobid_Custom_print_gmt=$(sbatch --job-name=$name_Custom_print_gmt --output=$outfile_Custom_print_gmt --partition=cpuq --time=24:00:00 --nodes=1 --ntasks-per-node=2 --mem-per-cpu=1024M --parsable --wrap="Rscript $Rscript_Custom_print_gmt --Table_of_gene_sets $Table_of_gene_sets --type $type --out $output_dir")
myjobid_seff_Custom_print_gmt=$(sbatch --dependency=afterany:$myjobid_Custom_print_gmt --open-mode=append --output=$outfile_Custom_print_gmt --job-name=$seff_name --partition=cpuq --time=24:00:00 --nodes=1 --ntasks-per-node=1 --mem-per-cpu=128M --parsable --wrap="seff $myjobid_Custom_print_gmt >> $outfile_Custom_print_gmt")

