#!/bin/bash

eval "$(conda shell.bash hook)"
  

Rscripts_path=$(echo "/home/manuel.tardaguila/Scripts/R/")

MASTER_ROUTE=$1
analysis=$2
annotation_file=$3
selected_annotations=$(echo "Dorothea_ABCD_GATA6_targets,Dorothea_ABCD_BCL11A_targets")
#selected_annotations_other=$(echo "CONCANNON_APOPTOSIS_BY_EPOXOMICIN_UP,G2M,HP_INCREASED_MEAN_PLATELET_VOLUME,HP_ABNORMAL_PLATELET_VOLUME")
selected_annotations_other=$(echo "ZHENG_CORD_BLOOD_C1_PUTATIVE_MEGAKARYOCYTE_PROGENITOR")


output_dir=$(echo "$MASTER_ROUTE""$analysis""/")

Log_files=$(echo "$output_dir""/""Log_files/")

conda activate multiome_NEW_downstream_analysis

type=$(echo "$analysis""_""bespoke_heatmap_function")
outfile_bespoke_heatmap_function=$(echo "$Log_files""outfile_6_""$type"".out")
touch $outfile_bespoke_heatmap_function
echo -n "" > $outfile_bespoke_heatmap_function
name_bespoke_heatmap_function=$(echo "$type""_job")
seff_name=$(echo "seff""_""$type")

Rscript_bespoke_heatmap_function=$(echo "$Rscripts_path""439_bespoke_DE_heatmap.R")

DE_result=$(echo "$output_dir""DE_results.rds")
normalised_counts=$(echo "$output_dir""norcounts_FINAL.rds")

## annotation_file=$(echo "$output_dir""genes_GSEA_annotated.tsv")


myjobid_bespoke_heatmap_function=$(sbatch --job-name=$name_bespoke_heatmap_function --output=$outfile_bespoke_heatmap_function --partition=cpuq --time=24:00:00 --nodes=1 --ntasks-per-node=2 --mem-per-cpu=4096 --parsable --wrap="Rscript $Rscript_bespoke_heatmap_function --DE_result $DE_result --normalised_counts $normalised_counts --selected_annotations $selected_annotations --annotation_file $annotation_file --selected_annotations_other $selected_annotations_other --type $type --out $output_dir")
myjobid_seff_bespoke_heatmap_function=$(sbatch --dependency=afterany:$myjobid_bespoke_heatmap_function --open-mode=append --output=$outfile_MSigDB_ORA --job-name=$seff_name --partition=cpuq --time=24:00:00 --nodes=1 --ntasks-per-node=1 --mem-per-cpu=128M --parsable --wrap="seff $myjobid_bespoke_heatmap_function >> $outfile_bespoke_heatmap_function")

conda deactivate
