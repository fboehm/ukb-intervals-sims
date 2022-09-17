#!/bin/bash

#SBATCH --partition=mulan,main
#SBATCH --time=1-00:00:00
#SBATCH --job-name=subset-plink-ref
#SBATCH --mem=2G
#SBATCH --output=/net/mulan/home/fredboe/research/ukb-intervals/cluster_outputs/subset-plink-ref.out
#SBATCH --error=/net/mulan/home/fredboe/research/ukb-intervals/cluster_outputs/subset-plink-ref.err

idxSub=~/research/ukb-intervals-sims/hapmap3/subjects_for_sims_reference.txt
let chr=1
      bfileAll=~/research/ukb-intervals-sims/hapmap3/plink_files_for_sims/chr${chr}
      bfileSub=~/research/ukb-intervals-sims/dat-quant/reference/allchr${chr}
      bfileSubP=~/research/ukb-intervals-sims/dat-quant/reference/chr${chr}
      plink-1.9 --silent --bfile ${bfileAll} --keep ${idxSub} --make-bed --out ${bfileSub}
      plink-1.9 --silent --bfile ${bfileSub} --maf 0.01 --make-bed --out ${bfileSubP}
      rm ${bfileSub}.bed
      rm ${bfileSub}.bim
      rm ${bfileSub}.fam
