#!/bin/bash

#SBATCH --partition=mulan,main
#SBATCH --time=1-00:00:00
#SBATCH --job-name=subset-plink-validation
#SBATCH --mem=2G
#SBATCH --output=/net/mulan/home/fredboe/research/ukb-intervals/cluster_outputs/subset-plink-val_%a.out
#SBATCH --error=/net/mulan/home/fredboe/research/ukb-intervals/cluster_outputs/subset-plink-val_%a.err


idxSub=~/research/ukb-intervals-sims/hapmap3-136k/subjects_for_sims_validation-136k.txt
let chr=1
      bfileAll=~/research/ukb-intervals-sims/hapmap3-136k/plink_files_for_sims/chr${chr}
      bfileSub=~/research/ukb-intervals-sims/dat-quant-136k/validation/allchr${chr}
      bfileSubP=~/research/ukb-intervals-sims/dat-quant-136k/validation/chr${chr}
      plink-1.9 --silent --bfile ${bfileAll} --keep ${idxSub} --make-bed --out ${bfileSub}
      plink-1.9 --silent --bfile ${bfileSub} --maf 0.01 --make-bed --out ${bfileSubP}
      rm ${bfileSub}.bed
      rm ${bfileSub}.bim
      rm ${bfileSub}.fam
