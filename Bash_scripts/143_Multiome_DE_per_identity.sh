#!/bin/bash

eval "$(conda shell.bash hook)"
  

Rscripts_path=$(echo "/home/manuel.tardaguila/Scripts/R/")

MASTER_ROUTE=$1
analysis=$2


output_dir=$(echo "$MASTER_ROUTE""$analysis""/")

Log_files=$(echo "$output_dir""/""Log_files/")

rm -rf $Log_files
mkdir -p $Log_files

conda activate multiome_NEW_downstream_analysis

### DE_function

type=$(echo "$analysis""_""DE_function")
outfile_DE_function=$(echo "$Log_files""outfile_1_""$type"".out")
touch $outfile_DE_function
echo -n "" > $outfile_DE_function
name_DE_function=$(echo "$type""_job")
seff_name=$(echo "seff""_""$type")

Rscript_DE_function=$(echo "$Rscripts_path""433_Multiome_DE_per_cell_type.R")

SeuratObject=$(echo "/group/soranzo/manuel.tardaguila/2025_hESC_MK_multiome/Downstream_analysis_cluster_after_genotyping/graphs_rpca/graphs_rpca/merged_clusters_after_genotyping_after_refined_annotation_new_peaks_rpca_integrate_rpca_annotation.rds")


myjobid_DE_function=$(sbatch --job-name=$name_DE_function --output=$outfile_DE_function --partition=cpuq --time=24:00:00 --nodes=1 --ntasks-per-node=8 --mem-per-cpu=4096 --parsable --wrap="Rscript $Rscript_DE_function --SeuratObject $SeuratObject --type $type --out $output_dir")
myjobid_seff_DE_function=$(sbatch --dependency=afterany:$myjobid_DE_function --open-mode=append --output=$outfile_DE_function --job-name=$seff_name --partition=cpuq --time=24:00:00 --nodes=1 --ntasks-per-node=1 --mem-per-cpu=128M --parsable --wrap="seff $myjobid_DE_function >> $outfile_DE_function")

### volcano_function

type=$(echo "$analysis""_""volcano_function")
outfile_volcano_function=$(echo "$Log_files""outfile_2_""$type"".out")
touch $outfile_volcano_function
echo -n "" > $outfile_volcano_function
name_volcano_function=$(echo "$type""_job")
seff_name=$(echo "seff""_""$type")

Rscript_volcano_function=$(echo "$Rscripts_path""434_Multiome_DE_volcano_plots.R")

DE_result=$(echo "$output_dir""DE_results.rds")

# --dependency=afterany:$myjobid_DE_function

myjobid_volcano_function=$(sbatch --dependency=afterany:$myjobid_DE_function --job-name=$name_volcano_function --output=$outfile_volcano_function --partition=cpuq --time=24:00:00 --nodes=1 --ntasks-per-node=2 --mem-per-cpu=4096 --parsable --wrap="Rscript $Rscript_volcano_function --DE_result $DE_result --type $type --out $output_dir")
myjobid_seff_volcano_function=$(sbatch --dependency=afterany:$myjobid_volcano_function --open-mode=append --output=$outfile_volcano_function --job-name=$seff_name --partition=cpuq --time=24:00:00 --nodes=1 --ntasks-per-node=1 --mem-per-cpu=128M --parsable --wrap="seff $myjobid_volcano_function >> $outfile_volcano_function")

### heatmap_function

type=$(echo "$analysis""_""heatmap_function")
outfile_heatmap_function=$(echo "$Log_files""outfile_3_""$type"".out")
touch $outfile_heatmap_function
echo -n "" > $outfile_heatmap_function
name_heatmap_function=$(echo "$type""_job")
seff_name=$(echo "seff""_""$type")

Rscript_heatmap_function=$(echo "$Rscripts_path""435_Multiome_DE_heatmap.R")

DE_result=$(echo "$output_dir""DE_results.rds")
normalised_counts=$(echo "$output_dir""norcounts_FINAL.rds")


# --dependency=afterany:$myjobid_DE_function

myjobid_heatmap_function=$(sbatch --dependency=afterany:$myjobid_DE_function --job-name=$name_heatmap_function --output=$outfile_heatmap_function --partition=cpuq --time=24:00:00 --nodes=1 --ntasks-per-node=2 --mem-per-cpu=4096 --parsable --wrap="Rscript $Rscript_heatmap_function --DE_result $DE_result --normalised_counts $normalised_counts --type $type --out $output_dir")
myjobid_seff_heatmap_function=$(sbatch --dependency=afterany:$myjobid_heatmap_function --open-mode=append --output=$outfile_heatmap_function --job-name=$seff_name --partition=cpuq --time=24:00:00 --nodes=1 --ntasks-per-node=1 --mem-per-cpu=128M --parsable --wrap="seff $myjobid_heatmap_function >> $outfile_heatmap_function")


conda deactivate
