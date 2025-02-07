#!/bin/bash
 
eval "$(conda shell.bash hook)"

output_dir=$1

Log_files=$(echo "$output_dir""/""Log_files/")

#rm -rf $Log_files
#mkdir -p $Log_files


declare -a arr

sample_array=$(echo "MCO_01326")

a=($(echo "$sample_array" | tr "," '\n'))

declare -a arr

for i  in "${a[@]}"
    do
        sample_array_sel=${i}
        echo "$sample_array_sel"


	module load cellranger/6.1.2


	type=$(echo "Run_cellranger""_""$sample_array_sel")


	outfile_Run_cellranger=$(echo "$Log_files""outfile_1_""$type"".log")
	touch $outfile_Run_cellranger
	echo -n "" > $outfile_Run_cellranger
	name_Run_cellranger=$(echo "$type""_job")





	#SBATCH --job-name=$name_Run_cellranger
	#SBATCH --mail-type=ALL
	#SBATCH --mail-user=manuel.tardaguila@fht.org
	#SBATCH --partition=cpuq
	#SBATCH --nodes=1
	#SBATCH --ntasks-per-node=16
	#SBATCH --output=$outfile_Run_cellranger
	#SBATCH --mem=48G
	#SBATCH --time=36:00:00


	cd $output_dir
	cellranger count --id $sample_array_sel \
		       --sample $sample_array_sel \
		       --transcriptome $(echo "/group/soranzo/manuel.tardaguila/2025_hESC_MK_multiome/GEX_reseq/alignment/GFP_transgene_vCHEK2_and_DNMT3A_cellranger") \
		       --fastqs $(echo "/group/soranzo/manuel.tardaguila/2025_hESC_MK_multiome/GEX_reseq/alignment/cellranger/") \
      		       --chemistry ARC-v1 \
		       --localcores=16 \
		       --localmem=47 \
		       --jobmode=local

	
done

