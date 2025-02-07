#!/bin/bash
 
eval "$(conda shell.bash hook)"
  

Rscripts_path=$(echo "/home/manuel.tardaguila/Scripts/R/")

MASTER_ROUTE=$1
analysis=$2

output_dir=$(echo "$MASTER_ROUTE""$analysis""/")



Log_files=$(echo "$output_dir""/""Log_files/")

#rm -rf $Log_files
#mkdir -p $Log_files

 

conda activate multiome_NEW_downstream_analysis

######################################################## LOOP #####################################################

declare -a arr

# problem early\ erythroid,early MK
cell_type_array=$(echo '1,2,3,4,5,6')

 a=($(echo "$cell_type_array" | tr "," '\n'))

 declare -a arr

 for i  in "${a[@]}"
 do

     cell_type_array_sel=$i
     echo "$cell_type_array_sel"

     ### DA_test

     type=$(echo "$cell_type_array_sel""_""DA_test")
     outfile_DA_test=$(echo "$Log_files""outfile_3_""$type"".log")
     touch $outfile_DA_test
     echo -n "" > $outfile_DA_test
     name_DA_test=$(echo "$type""_job")
     seff_name=$(echo "seff""_""$type")



     Rscript_DA_test=$(echo "$Rscripts_path""418_Per_cell_type_PB_ATAC_CHEK2.R")

     adata=$(echo "/group/soranzo/manuel.tardaguila/2025_hESC_MK_multiome/processing_outputs/merged_processed_presentation_updated_integrative_genotyping.rds")
     cell_type_sel=$cell_type_array_sel
     
     myjobid_DA_test=$(sbatch --job-name $name_DA_test --output=$outfile_DA_test --partition=cpuq --time=24:00:00 --nodes=1 --ntasks-per-node=8 --mem-per-cpu=8096 --parsable --wrap="Rscript $Rscript_DA_test --adata $adata --cell_type_sel $cell_type_sel --type $type --out $output_dir")
     myjobid_seff_DA_test=$(sbatch --dependency=afterany:$myjobid_DA_test --open-mode=append --output=$outfile_DA_test --job-name=$seff_name --partition=cpuq --time=24:00:00 --nodes=1 --ntasks-per-node=1 --mem-per-cpu=128M --parsable --wrap="seff $myjobid_DA_test >> $outfile_DA_test")

    
     echo "->>>$myjobid_DA_test"
     arr[${#arr[@]}]="$myjobid_DA_test"

     
 done
      


conda deactivate
