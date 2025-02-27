#!/bin/bash

eval "$(conda shell.bash hook)"
  

Rscripts_path=$(echo "/home/manuel.tardaguila/Scripts/R/")

MASTER_ROUTE=$1
analysis=$2
DE_results=$3
TF_terms=$(echo 'RUNX1,GATA2,CUX1,CREB5,GATA1,GATA6,BCL11A,BCL11B,MAZ,NR2F1,PROX1,ZBTB7A')
search_terms=$(echo 'PLATELET,ERYTHROCYTE,MEGAKARYOCYTE,CHEK2,ATM,DNMT3A,HEMATOPOIETIC,HEMATOPOIESIS,CLONAL,EMBRYONIC,STEM,GRAHAM_NORMAL_QUIESCENT_VS_NORMAL_DIVIDING_DN,RUNX1,GATA2,CUX1,CREB5,GATA1,GATA6,BCL11A,BCL11B,MAZ,NR2F1,PROX1,ZBTB7A')
path_to_GMT=$(echo "/home/manuel.tardaguila/GMT_files/msigdb_v2023.1.Hs_files_to_download_locally_ENTREZ/")

# DE_results=$(echo "/group/soranzo/manuel.tardaguila/2025_hESC_MK_multiome/Downstream_analysis/DE_results_without_time_as_a_covariate_NEW_METHOD.tsv")

output_dir=$(echo "$MASTER_ROUTE""$analysis""/")

Log_files=$(echo "$output_dir""/""Log_files/")

rm -rf $Log_files
mkdir -p $Log_files


### MSigDB_ORA


type=$(echo "$analysis""_""MSigDB_ORA")
outfile_MSigDB_ORA=$(echo "$Log_files""outfile_1_""$type"".out")
touch $outfile_MSigDB_ORA
echo -n "" > $outfile_MSigDB_ORA
name_MSigDB_ORA=$(echo "$type""_job")
seff_name=$(echo "seff""_""$type")

#Rscript_MSigDB_ORA=$(echo "$Rscripts_path""419_MySigDb_DESeq2_on_MK_diff_v_cluster.R")
Rscript_MSigDB_ORA=$(echo "$Rscripts_path""419_MySigDb_DESeq2_on_MK_diff_v_cell_type.R")


background_genes=$(echo "/home/manuel.tardaguila/GMT_files/msigdb_v2023.1.Hs_files_to_download_locally_ENTREZ/c5.hpo.v2024.1.Hs.entrez.gmt")
maxGS_size=$(echo "500")
minGS_size=$(echo "10")
pval_threshold=$(echo "0.05")
log2FC_threshold=$(echo "0.25")
seff_name=$(echo "seff""_""$type")seff_name=$(echo "seff""_""$type")

myjobid_MSigDB_ORA=$(sbatch --job-name=$name_MSigDB_ORA --output=$outfile_MSigDB_ORA --partition=cpuq --time=24:00:00 --nodes=1 --ntasks-per-node=2 --mem-per-cpu=1024 --parsable --wrap="Rscript $Rscript_MSigDB_ORA --DE_results $DE_results --path_to_GMT $path_to_GMT --search_terms $search_terms --background_genes $background_genes --exemption_terms $TF_terms --maxGS_size $maxGS_size --minGS_size $minGS_size --pval_threshold $pval_threshold --log2FC_threshold $log2FC_threshold --type $type --out $output_dir")
myjobid_seff_MSigDB_ORA=$(sbatch --dependency=afterany:$myjobid_MSigDB_ORA --open-mode=append --output=$outfile_MSigDB_ORA --job-name=$seff_name --partition=cpuq --time=24:00:00 --nodes=1 --ntasks-per-node=1 --mem-per-cpu=128M --parsable --wrap="seff $myjobid_MSigDB_ORA >> $outfile_MSigDB_ORA")

### MSigDB_GSEA

conda activate GSEA

type=$(echo "$analysis""_""MSigDB_GSEA")
outfile_MSigDB_GSEA=$(echo "$Log_files""outfile_2_""$type"".out")
touch $outfile_MSigDB_GSEA
echo -n "" > $outfile_MSigDB_GSEA
name_MSigDB_GSEA=$(echo "$type""_job")
seff_name=$(echo "seff""_""$type")

#Rscript_MSigDB_GSEA=$(echo "$Rscripts_path""420_GSEA_DESeq2_on_MK_diff_v_cluster.R")
Rscript_MSigDB_GSEA=$(echo "$Rscripts_path""420_GSEA_DESeq2_on_MK_diff_v_cell_type.R")



pval_threshold=$(echo "0.05")
log2FC_threshold=$(echo "0")
Threshold_number_of_genes=$(echo '3')


myjobid_MSigDB_GSEA=$(sbatch --job-name=$name_MSigDB_GSEA --output=$outfile_MSigDB_GSEA --partition=cpuq --time=24:00:00 --nodes=1 --ntasks-per-node=2 --mem-per-cpu=1024 --parsable --wrap="Rscript $Rscript_MSigDB_GSEA --DE_results $DE_results --path_to_GMT $path_to_GMT --search_terms $search_terms --pval_threshold $pval_threshold --log2FC_threshold $log2FC_threshold --Threshold_number_of_genes $Threshold_number_of_genes --TF_terms $TF_terms --type $type --out $output_dir")
myjobid_seff_MSigDB_GSEA=$(sbatch --dependency=afterany:$myjobid_MSigDB_GSEA --open-mode=append --output=$outfile_MSigDB_GSEA --job-name=$seff_name --partition=cpuq --time=24:00:00 --nodes=1 --ntasks-per-node=1 --mem-per-cpu=128M --parsable --wrap="seff $myjobid_MSigDB_GSEA >> $outfile_MSigDB_GSEA")
 


conda deactivate
