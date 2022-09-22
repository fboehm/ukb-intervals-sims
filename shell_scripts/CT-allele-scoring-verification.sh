#!/bin/bash

#SBATCH --partition=main
#SBATCH --time=1:00:00
#SBATCH --job-name=allele-verif
#SBATCH --mem=2G
#SBATCH --array=1-450%100
#SBATCH --output=/net/mulan/home/fredboe/research/ukb-intervals-sims/cluster_outputs/allele-scoring-verif_%a.out
#SBATCH --error=/net/mulan/home/fredboe/research/ukb-intervals-sims/cluster_outputs/allele-scoring-verif_%a.err

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
idxverif=~/research/ukb-intervals-sims/hapmap3/subjects_for_sims_verification.txt
#esteffdbslmmt=${compstr}05_internal_c/pheno${p}/DBSLMM/summary_hm3_cross${cross}_chr${chr}_best.dbslmm.txt
esteffdbslmmt=~/research/ukb-intervals-sims/dat-quant/CT/scenario${scenario}/${distribution}/hsq${hsq}/esteff_ukb_pheno${p}_cross${fold}.txt
outpath=~/research/ukb-intervals-sims/dat-quant/verification/CT-allele-scores/scenario${scenario}/${distribution}/hsq${hsq}
mkdir -p ${outpath}
preddbslmmt=${outpath}/pred_replicate${p}_fold${fold}.txt
outfile=${preddbslmmt}.profile
#gunzip ${esteffdbslmmt}.gz
if [ ! -f "$outfile" ]; then # check if file doesn't exist
plink-1.9 --silent --bfile ${bfile} --score ${esteffdbslmmt} 1 2 3 sum --keep ${idxverif} --out ${preddbslmmt} --allow-no-sex
# plink-1.9 --silent --bfile ${bfile} --score ${esteffdbslmmt} 1 2 4 sum --keep ${idxagg} --out ${aggdbslmm
fi
fi
done
done
done
done
done
