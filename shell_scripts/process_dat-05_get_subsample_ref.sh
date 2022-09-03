#!/bin/bash


#SBATCH --partition=mulan,main
#SBATCH --time=1-00:00:00
#SBATCH --job-name=subset-ref
#SBATCH --mem-per-cpu=2G
#SBATCH --array=1-22
#SBATCH --output=/net/mulan/home/fredboe/research/ukb-intervals-sims/cluster_outputs/process_dat-05_get_subsample_ref_%j_%a.out
#SBATCH --error=/net/mulan/home/fredboe/research/ukb-intervals-sims/cluster_outputs/process_dat-05_get_subsample_ref_%j_%a.err


compStr=/net/mulan/disk2/yasheng/comparisonProject/
fbStr=~/research/ukb-intervals-sims/dat/
# ukb subsample
let k=0
for chr in `seq 1 22`;do
let k=${k}+1
if [ ${k} -eq ${SLURM_ARRAY_TASK_ID} ]; then


# ## subsample
bfileAll=/net/mulan/disk2/yasheng/predictionProject/plink_file/ukb/chr${chr}
## reference sample
bfileRef=${fbStr}reference/ukb/geno/allchr${chr}
bfileRefP=${fbStr}reference/ukb/geno/chr${chr}
idxRef=${fbStr}reference/01_idx.txt
plink-1.9 --bfile ${bfileAll} --keep ${idxRef} --make-bed --out ${bfileRef}
plink-1.9 --bfile ${bfileRef} --maf 0.01 --make-bed --out ${bfileRefP}

## remove files
rm ${bfileRef}.log
rm ${bfileRef}.bed
rm ${bfileRef}.bim
rm ${bfileRef}.fam
rm ${bfileRefP}.log

done


fi