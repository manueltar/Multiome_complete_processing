#!/bin/bash

eval "$(conda shell.bash hook)"
 
  
Rscripts_path=$(echo "/home/manuel.tardaguila/Scripts/R/simba/")

MASTER_ROUTE=$1
analysis=$2


##########################################################################

output_dir=$(echo "$MASTER_ROUTE""$analysis""/")

Log_files=$(echo "$output_dir""/""Log_files/")
 
conda activate SIMBA_R
module load bedtools2/2.31.0
 
### simba_scan_for_kmers_motifs

type=$(echo "simba_scan_for_kmers_motifs")
outfile_simba_scan_for_kmers_motifs=$(echo "$Log_files""outfile_8_""$type"".log")
touch $outfile_simba_scan_for_kmers_motifs
echo -n "" > $outfile_simba_scan_for_kmers_motifs
name_simba_scan_for_kmers_motifs=$(echo "$type""_job")


Rscript_simba_scan_for_kmers_motifs=$(echo "$Rscripts_path""scan_for_kmers_motifs.R")


input_bed=$(echo "/group/soranzo/manuel.tardaguila/2025_hESC_MK_multiome/Downstream_analysis_cluster_after_genotyping/Peaks.bed")
reference_genome=$(echo "/processing_data/reference_datasets/iGenomes/2022.1/Homo_sapiens/NCBI/GRCh38/Sequence/WholeGenomeFasta/genome.fa")
species=$(echo 'Homo_sapiens')
outdir=$(echo "$output_dir""result_SIMBA")

myjobid_simba_scan_for_kmers_motifs=$(sbatch --job-name $name_simba_scan_for_kmers_motifs --output=$outfile_simba_scan_for_kmers_motifs --partition=cpuq --time=24:00:00 --nodes=1 --ntasks-per-node=20 --mem-per-cpu=4096 --parsable --wrap="Rscript $Rscript_simba_scan_for_kmers_motifs --input $input_bed --genome $reference_genome --species $species --output $outdir")
myjobid_seff_simba_scan_for_kmers_motifs=$(sbatch --dependency=afterany:$myjobid_simba_scan_for_kmers_motifs --open-mode=append --output=$outfile_simba_scan_for_kmers_motifs --job-name="seff" --partition=cpuq --time=24:00:00 --nodes=1 --ntasks-per-node=1 --mem-per-cpu=128M --parsable --wrap="seff $myjobid_simba_scan_for_kmers_motifs >> $outfile_simba_scan_for_kmers_motifs")

conda deactivate
