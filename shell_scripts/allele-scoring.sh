#!/bin/bash

#SBATCH --partition=mulan,main
#SBATCH --time=1-00:00:00
#SBATCH --job-name=allele-scoring
#SBATCH --mem=2G
#SBATCH --array=1-1800%100
#SBATCH --output=/net/mulan/home/fredboe/research/ukb-intervals-sims/cluster_outputs/allele-scoring_%a.out
#SBATCH --error=/net/mulan/home/fredboe/research/ukb-intervals-sims/cluster_outputs/allele-scoring_%a.err


scenarios=( I II III IV)
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

esteffdbslmmt=~/research/ukb-intervals-sims/dat-quant/DBSLMM/scenario${scenario}/${distribution}/hsq${hsq}/summary_ukb_pheno${p}_scenario${scenario}_${distribution}_hsq${hsq}_fold${fold}_chr${chr}_best.dbslmm.txt
preddbslmmt=~/research/ukb-intervals-sims/dat-quant/DBSLMM/scenario${scenario}/${distribution}/hsq${hsq}/pred_ukb_pheno${p}_scenario${scenario}_${distribution}_hsq${hsq}_fold${fold}_chr${chr}_best.dbslmm.txt
plink-1.9 --silent --bfile ${bfile} --score ${esteffdbslmmt} 1 2 4 sum --keep ${idxtest} --out ${preddbslmmt} --allow-no-sex
fi
done
done
done
done
done