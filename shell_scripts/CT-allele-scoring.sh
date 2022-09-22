#!/bin/bash

#SBATCH --partition=mulan,main
#SBATCH --time=1-00:00:00
#SBATCH --job-name=CT-allele-scoring
#SBATCH --mem=2G
#SBATCH --array=1-450%100
#SBATCH --output=/net/mulan/home/fredboe/research/ukb-intervals-sims/cluster_outputs/CT-allele-scoring_%a.out
#SBATCH --error=/net/mulan/home/fredboe/research/ukb-intervals-sims/cluster_outputs/CT-allele-scoring_%a.err


scenarios=( I )
distributions=( laplace normal scaledt)
hsqs=(0.1 0.2 0.5)


# 
let k=0
#let pc_ctr=0
chr=1


for scenario in ${scenarios[@]}; do
  for distribution in ${distributions[@]}; do
    for hsq in ${hsqs[@]}; do
      for p in `seq 1 10`; do
        for fold in `seq 1 5`; do

let k=${k}+1
if [ ${k} -eq ${SLURM_ARRAY_TASK_ID} ]
then
# bfile
bfile=~/research/ukb-intervals-sims/hapmap3/plink_files_for_sims/chr${chr}
idxtest=~/research/ukb-intervals-sims/hapmap3/test-ids-fold${fold}.txt

esteffdbslmmt=~/research/ukb-intervals-sims/dat-quant/CT/scenario${scenario}/${distribution}/hsq${hsq}/esteff_ukb_pheno${p}_cross${fold}.txt
preddbslmmt=~/research/ukb-intervals-sims/dat-quant/CT/scenario${scenario}/${distribution}/hsq${hsq}/pred_ukb_pheno${p}_fold${fold}.txt
outfile=${preddbslmmt}.profile
#if [ ! -f "$outfile" ]; then # check if file doesn't exist
plink-1.9 --silent --bfile ${bfile} --score ${esteffdbslmmt} 1 2 3 sum --keep ${idxtest} --out ${preddbslmmt} --allow-no-sex
#fi
fi
done
done
done
done
done