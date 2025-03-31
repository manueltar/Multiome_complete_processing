#!/usr/bin/env bash
#
# =============================================================================
# Job Script
# =============================================================================
#

#SBATCH --job-name=merge_frag
#SBATCH --mail-type=ALL
#SBATCH --mail-user=manuel.tardaguila@fht.org
#SBATCH --partition=cpuq
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=16
#SBATCH --output=merge_frag_%j.log
#SBATCH --mem=96G

# 16 96 --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------


echo "Running on $(hostname)"
echo "Started: $(date)"
echo "========================"

MASTER_ROUTE=$1

module load htslib-tools

cd $MASTER_ROUTE

for SAMPLE in MCO_01326 MCO_01327 MCO_01328 MCO_01329 MCO_01330;
do echo $SAMPLE

   
   if [[ ($SAMPLE = 'MCO_01326') || ($SAMPLE = 'MCO_01327') || ($SAMPLE = 'MCO_01328') || ($SAMPLE = 'MCO_01329') ]]; then

       FILE=$(echo "/group/soranzo/manuel.tardaguila/2025_hESC_MK_multiome/""$SAMPLE""/outs/atac_fragments.tsv.gz")



   else
       
       FILE=$(echo "/group/soranzo/manuel.tardaguila/2025_hESC_lymph_multiome/Multiome/""$SAMPLE""/outs/atac_fragments.tsv.gz")
   fi

   echo "$FILE"

   zcat $FILE | awk -v SAMPLE=$SAMPLE \ 'BEGIN{FS=OFS="\t"} {print $1,$2,$3,SAMPLE"_"$4,$5}\';
done | sort -k1,1 -k2,2n -S 72G | bgzip -c -@ 16 > $MASTER_ROUTE/processing_outputs/merged.atac_fragments.tsv.gz 

tabix -p bed $MASTER_ROUTE/processing_outputs/merged.atac_fragments.tsv.gz



echo "========================"
echo "Completed: $(date)"
