#!/bin/bash

eval "$(conda shell.bash hook)"
 

python_path=$(echo "/home/manuel.tardaguila/Scripts/Python/")

MASTER_ROUTE=$1
analysis=$2

output_dir=$(echo "$MASTER_ROUTE""$analysis""/")

Log_files=$(echo "$output_dir""/""Log_files/")

echo "$Log_files"

rm -rf $Log_files
mkdir -p $Log_files

conda activate /group/soranzo/conda_envs/env_simba

### SIMBA_QC_and_embeddings

type=$(echo "SIMBA_QC_and_embeddings")
outfile_SIMBA_QC_and_embeddings=$(echo "$Log_files""outfile_1_""$type"".log")
touch $outfile_SIMBA_QC_and_embeddings
echo -n "" > $outfile_SIMBA_QC_and_embeddings
name_SIMBA_QC_and_embeddings=$(echo "$type""_job")


python_script_SIMBA_QC_and_embeddings=$(echo "$python_path""1_simba_QC_model_training_all_features.py")

memory=$(echo "12000")
processors=$(echo "30")

#memory=$(echo "8192")
#processors=$(echo "10")

### I had to edit the file /group/soranzo/manuel.tardaguila/2025_hESC_MK_multiome/Downstream_analysis_cluster_after_genotyping/result_SIMBA/pbg/graph0/model/config.json to include the full path 

cd $MASTER_ROUTE

input_rna=$(echo "/group/soranzo/manuel.tardaguila/2025_hESC_MK_multiome/Downstream_analysis_cluster_after_genotyping/RNA.h5ad")
input_atac=$(echo "/group/soranzo/manuel.tardaguila/2025_hESC_MK_multiome/Downstream_analysis_cluster_after_genotyping/ATAC.h5ad")
input_kmers=$(echo "/group/soranzo/manuel.tardaguila/2025_hESC_MK_multiome/Downstream_analysis_cluster_after_genotyping/result_SIMBA/freq_kmer.h5")
input_motifs=$(echo "/group/soranzo/manuel.tardaguila/2025_hESC_MK_multiome/Downstream_analysis_cluster_after_genotyping/result_SIMBA//freq_motif.h5")
preexisting_graph=$(echo "/group/soranzo/manuel.tardaguila/2025_hESC_MK_multiome/Downstream_analysis_cluster_after_genotyping/result_SIMBA/pbg/graph0/")
preexisting_model=$(echo "/group/soranzo/manuel.tardaguila/2025_hESC_MK_multiome/Downstream_analysis_cluster_after_genotyping/result_SIMBA/pbg/graph0/model/")
path_entity_alias=$(echo "/group/soranzo/manuel.tardaguila/2025_hESC_MK_multiome/Downstream_analysis_cluster_after_genotyping/result_SIMBA/pbg/graph0/")
path_entity=$(echo "/group/soranzo/manuel.tardaguila/2025_hESC_MK_multiome/Downstream_analysis_cluster_after_genotyping/result_SIMBA/pbg/graph0/input/entity/")
path_emb=$(echo "/group/soranzo/manuel.tardaguila/2025_hESC_MK_multiome/Downstream_analysis_cluster_after_genotyping/result_SIMBA/pbg/graph0/model/")

myjobid_SIMBA_QC_and_embeddings=$(sbatch --job-name $name_SIMBA_QC_and_embeddings --output=$outfile_SIMBA_QC_and_embeddings --partition=cpuq --time=24:00:00 --nodes=1 --ntasks-per-node=$processors --mem-per-cpu=$memory --parsable --wrap="python3 $python_script_SIMBA_QC_and_embeddings --threads $processors --input-rna $input_rna --input-atac $input_atac --input-kmers $input_kmers --input-motifs $input_motifs --skip-graph-generation --skip-train-model --preexisting-graph $preexisting_graph --preexisting-model $preexisting_model --path-entity_alias $path_entity_alias --path-entity $path_entity --path-emb $path_emb --output $output_dir")
myjobid_seff_SIMBA_QC_and_embeddings=$(sbatch --dependency=afterany:$myjobid_SIMBA_QC_and_embeddings --open-mode=append --output=$outfile_SIMBA_QC_and_embeddings --job-name="seff" --partition=cpuq --time=24:00:00 --nodes=1 --ntasks-per-node=1 --mem-per-cpu=128M --parsable --wrap="seff $myjobid_SIMBA_QC_and_embeddings >> $outfile_SIMBA_QC_and_embeddings")

conda deactivate


# --skip-graph-generation --skip-train-model --skip-kmers_motifs --skip-compare_entities --skip-global_embedding --skip-RNA_QC --skip-ATAC_QC --skip-global_embedding --skip-compare_entities
