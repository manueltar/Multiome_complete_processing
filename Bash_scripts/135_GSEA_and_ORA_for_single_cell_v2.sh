#!/bin/bash

eval "$(conda shell.bash hook)"
  

Rscripts_path=$(echo "/home/manuel.tardaguila/Scripts/R/")

MASTER_ROUTE=$1
analysis=$2
DE_results=$3
TF_terms=$(echo 'RUNX1,GATA2,CUX1,CREB5,GATA1,GATA6,BCL11A,BCL11B,MAZ,NR2F1,PROX1,ZBTB7A')
search_terms=$(echo 'FIBROSIS,PLATELET,ERYTHROCYTE,MEGAKARYOCYTE,CHEK2,_ATM_,DNMT3A,HEMATOPOIETIC_,HEMATOPOIESIS_,RUNX1,GATA2,CUX1,CREB5,GATA1,GATA6,BCL11A,BCL11B,MAZ,NR2F1,PROX1,ZBTB7A,GFI,GFI1B')
path_to_GMT=$(echo "/home/manuel.tardaguila/GMT_files/msigdb_v2023.1.Hs_files_to_download_locally_ENTREZ/")



output_dir=$(echo "$MASTER_ROUTE""$analysis""/")

Log_files=$(echo "$output_dir""/""Log_files/")

rm -rf $Log_files
mkdir -p $Log_files

conda activate GSEA

### MSigDB_ORA

type=$(echo "$analysis""_""MSigDB_ORA")
outfile_MSigDB_ORA=$(echo "$Log_files""outfile_1_""$type"".out")
touch $outfile_MSigDB_ORA
echo -n "" > $outfile_MSigDB_ORA
name_MSigDB_ORA=$(echo "$type""_job")
seff_name=$(echo "seff""_""$type")

Rscript_MSigDB_ORA=$(echo "$Rscripts_path""419_ORA_DESeq2_on_MK_diff_v_cell_type_vclusterProfiler.R")



pval_threshold=$(echo "0.05")
log2FC_threshold=$(echo "0")
Threshold_number_of_genes=$(echo '3')


myjobid_MSigDB_ORA=$(sbatch --job-name=$name_MSigDB_ORA --output=$outfile_MSigDB_ORA --partition=cpuq --time=24:00:00 --nodes=1 --ntasks-per-node=2 --mem-per-cpu=1024 --parsable --wrap="Rscript $Rscript_MSigDB_ORA --DE_results $DE_results --path_to_GMT $path_to_GMT --search_terms $search_terms --pval_threshold $pval_threshold --log2FC_threshold $log2FC_threshold --Threshold_number_of_genes $Threshold_number_of_genes --TF_terms $TF_terms --type $type --out $output_dir")
myjobid_seff_MSigDB_ORA=$(sbatch --dependency=afterany:$myjobid_MSigDB_ORA --open-mode=append --output=$outfile_MSigDB_ORA --job-name=$seff_name --partition=cpuq --time=24:00:00 --nodes=1 --ntasks-per-node=1 --mem-per-cpu=128M --parsable --wrap="seff $myjobid_MSigDB_ORA >> $outfile_MSigDB_ORA")

### MSigDB_GSEA

type=$(echo "$analysis""_""MSigDB_GSEA")
outfile_MSigDB_GSEA=$(echo "$Log_files""outfile_2_""$type"".out")
touch $outfile_MSigDB_GSEA
echo -n "" > $outfile_MSigDB_GSEA
name_MSigDB_GSEA=$(echo "$type""_job")
seff_name=$(echo "seff""_""$type")

Rscript_MSigDB_GSEA=$(echo "$Rscripts_path""420_GSEA_DESeq2_on_MK_diff_v_cell_type_vclusterProfiler.R")



pval_threshold=$(echo "0.05")
log2FC_threshold=$(echo "0")
Threshold_number_of_genes=$(echo '3')


myjobid_MSigDB_GSEA=$(sbatch --job-name=$name_MSigDB_GSEA --output=$outfile_MSigDB_GSEA --partition=cpuq --time=24:00:00 --nodes=1 --ntasks-per-node=2 --mem-per-cpu=1024 --parsable --wrap="Rscript $Rscript_MSigDB_GSEA --DE_results $DE_results --path_to_GMT $path_to_GMT --search_terms $search_terms --pval_threshold $pval_threshold --log2FC_threshold $log2FC_threshold --Threshold_number_of_genes $Threshold_number_of_genes --TF_terms $TF_terms --type $type --out $output_dir")
myjobid_seff_MSigDB_GSEA=$(sbatch --dependency=afterany:$myjobid_MSigDB_GSEA --open-mode=append --output=$outfile_MSigDB_GSEA --job-name=$seff_name --partition=cpuq --time=24:00:00 --nodes=1 --ntasks-per-node=1 --mem-per-cpu=128M --parsable --wrap="seff $myjobid_MSigDB_ORA >> $outfile_MSigDB_ORA")



conda deactivate
