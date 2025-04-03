#!/bin/bash

eval "$(conda shell.bash hook)"
 

Rscripts_path=$(echo "/home/manuel.tardaguila/Scripts/R/")

conda activate multiome_QC_DEF

#sample_array=$(echo 'MCO_01326,MCO_01327,MCO_01328,MCO_01329,MCO_01330')
sample_array=$(echo 'MCO_01330,MCO_01331,MCO_01332,MCO_01333')


a=($(echo "$sample_array" | tr "," '\n'))

declare -a arr

for i  in "${a[@]}"
do

    sample_array_sel=$i
    echo "$sample_array_sel"

##   rm -rf $premerge_dir ######################################################################### SECURITY
##   mkdir -p $premerge_dir ####################################################################### SECURITY

     if [[ ($sample_array_sel = 'MCO_01326') || ($sample_array_sel = 'MCO_01327') || ($sample_array_sel = 'MCO_01328') || ($sample_array_sel = 'MCO_01329') ]]; then

	 preliminary_filtered=$(echo "/group/soranzo/manuel.tardaguila/2025_hESC_MK_multiome/""processing_outputs""/""$sample_array_sel""/""intermediate""/""preliminary_filtered.rds")
	 path_processing_outputs=$(echo "/group/soranzo/manuel.tardaguila/2025_hESC_MK_multiome/""processing_outputs""/""$sample_array_sel""/")
	 intermediate_dir=$(echo "/group/soranzo/manuel.tardaguila/2025_hESC_MK_multiome/""processing_outputs""/""$sample_array_sel""/""intermediate""/")
	 snATAC_dir=$(echo "/group/soranzo/manuel.tardaguila/2025_hESC_MK_multiome/""processing_outputs""/""$sample_array_sel""/""snATAC_matrices""/")
	 crange_dir=$(echo "/group/soranzo/manuel.tardaguila/2025_hESC_MK_multiome/""$sample_array_sel""/""outs""/")
	 premerge_dir=$(echo "/group/soranzo/manuel.tardaguila/2025_hESC_MK_multiome/""processing_outputs""/""$sample_array_sel""/""pre_merge""/")
	 fragfile=$(echo "/group/soranzo/manuel.tardaguila/2025_hESC_MK_multiome/""$sample_array_sel""/""outs""/""atac_fragments.tsv.gz")

	 output_dir=$(echo "/group/soranzo/manuel.tardaguila/2025_hESC_MK_multiome/""processing_outputs""/")
	 Log_files=$(echo "$output_dir""/""Log_files/")

   else
       	 preliminary_filtered=$(echo "/group/soranzo/manuel.tardaguila/2025_hESC_lymph_multiome/Multiome/""processing_outputs""/""$sample_array_sel""/""intermediate""/""preliminary_filtered.rds")
	 path_processing_outputs=$(echo "/group/soranzo/manuel.tardaguila/2025_hESC_lymph_multiome/Multiome/""processing_outputs""/""$sample_array_sel""/")
	 intermediate_dir=$(echo "/group/soranzo/manuel.tardaguila/2025_hESC_lymph_multiome/Multiome/""processing_outputs""/""$sample_array_sel""/""intermediate""/")
	 snATAC_dir=$(echo "/group/soranzo/manuel.tardaguila/2025_hESC_lymph_multiome/Multiome/""processing_outputs""/""$sample_array_sel""/""snATAC_matrices""/")
	 crange_dir=$(echo "/group/soranzo/manuel.tardaguila/2025_hESC_lymph_multiome/Multiome/""$sample_array_sel""/""outs""/")
	 premerge_dir=$(echo "/group/soranzo/manuel.tardaguila/2025_hESC_lymph_multiome/Multiome/""processing_outputs""/""$sample_array_sel""/""pre_merge""/")
 	 fragfile=$(echo "/group/soranzo/manuel.tardaguila/2025_hESC_lymph_multiome/Multiome/""$sample_array_sel""/""outs""/""atac_fragments.tsv.gz")

	 output_dir=$(echo "/group/soranzo/manuel.tardaguila/2025_hESC_lymph_multiome/Multiome/""processing_outputs""/")
	 Log_files=$(echo "$output_dir""/""Log_files/")
	 
   fi

   echo "$preliminary_filtered"
   echo "$path_processing_outputs"
   echo "$intermediate_dir"
   echo "$snATAC_dir"
   echo "$crange_dir"
   echo "$premerge_dir"
   echo "$fragfile"

   
   rm -rf $premerge_dir
   mkdir -p $premerge_dir

    ### Seurat_second_pass

    type=$(echo "$sample_array_sel""_""Seurat_second_pass")
    outfile_Seurat_second_pass=$(echo "$Log_files""outfile_5_""$type"".log")
    touch $outfile_Seurat_second_pass
    echo -n "" > $outfile_Seurat_second_pass
    name_Seurat_second_pass=$(echo "$type""_job")

 
    Rscript_Seurat_second_pass=$(echo "$Rscripts_path""405_Seurat_second_pass_v3.R")

    sample_name=$sample_array_sel
 

    myjobid_Seurat_second_pass=$(sbatch --job-name $name_Seurat_second_pass --output=$outfile_Seurat_second_pass --partition=cpuq --time=24:00:00 --nodes=1 --ntasks-per-node=15 --mem-per-cpu=8192 --parsable --wrap="Rscript $Rscript_Seurat_second_pass --sample_name $sample_name --preliminary_filtered $preliminary_filtered --path_processing_outputs $path_processing_outputs --intermediate_dir $intermediate_dir --snATAC_dir $snATAC_dir --crange_dir $crange_dir --premerge_dir $premerge_dir --fragfile $fragfile  --type $type --out $output_dir")
    myjobid_seff_Seurat_second_pass=$(sbatch --dependency=afterany:$myjobid_Seurat_second_pass --open-mode=append --output=$outfile_Seurat_second_pass --job-name="seff" --partition=cpuq --time=24:00:00 --nodes=1 --ntasks-per-node=1 --mem-per-cpu=128M --parsable --wrap="seff $myjobid_Seurat_second_pass >> $outfile_Seurat_second_pass")


done

conda deactivate
