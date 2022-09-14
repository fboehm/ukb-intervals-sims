#!/bin/bash

#SBATCH --partition=mulan,main
#SBATCH --time=1-00:00:00
#SBATCH --job-name=plink_subset_SNPs
#SBATCH --mem=64G
#SBATCH --output=/net/mulan/home/fredboe/research/ukb-intervals-sims/cluster_outputs/subset_snps.out
#SBATCH --error=/net/mulan/home/fredboe/research/ukb-intervals-sims/cluster_outputs/subset_snps.err


dir=~/research/ukb-intervals-sims/hapmap3/ukb/
snplist=~/research/ukb-intervals-sims/hapmap3/snp_list.txt


#let k=0
for chr in `seq 1 21`; do
#  let k=${k}+1
#  if [ ${k} -eq ${SLURM_ARRAY_TASK_ID} ]; then
    plink --bfile ${dir}chr${chr} --extract ${snplist} --make-bed --allow-no-sex --out ~/research/ukb-intervals-sims/hapmap3/hm3_ukb/chr${chr}
#    fi
  done
