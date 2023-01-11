#!/bin/bash

#SBATCH --partition=mulan,main
#SBATCH --time=8:00:00
#SBATCH --job-name=plink_subset_SNPs
#SBATCH --mem=16G
#SBATCH --output=/net/mulan/home/fredboe/research/ukb-intervals-sims/cluster_outputs/subset_snps.out
#SBATCH --error=/net/mulan/home/fredboe/research/ukb-intervals-sims/cluster_outputs/subset_snps.err

# GOAL: Create a set of plink files for Chr1 for use in simulating traits
dir=~/research/ukb-intervals-sims/hapmap3/hm3_ukb/
subject_list=~/research/ukb-intervals-sims/hapmap3-small-scale-simulations/subjects_for_sims.txt
outdir=~/research/ukb-intervals-sims/hapmap3-small-scale-simulations/plink_files_for_sims/

mkdir -p ${outdir}

for chr in 1; do
    plink --bfile ${dir}chr${chr} --keep ${subject_list} --make-bed --allow-no-sex --out ${outdir}chr${chr}
done
