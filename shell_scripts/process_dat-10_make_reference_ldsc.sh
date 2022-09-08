#!/bin/bash

#SBATCH --partition=mulan,main
#SBATCH --time=1-00:00:00
#SBATCH --job-name=ldsc-ref
#SBATCH --mem-per-cpu=6G
#SBATCH --array=1-22
#SBATCH --output=/net/mulan/home/fredboe/research/ukb-intervals-sims/cluster_outputs/process_dat-10_%j_%a.out
#SBATCH --error=/net/mulan/home/fredboe/research/ukb-intervals-sims/cluster_outputs/process_dat-10_%j_%a.err


compStr=/net/mulan/disk2/yasheng/comparisonProject/
fbStr=~/research/ukb-intervals-sims/dat/

let k=0
for chr in `seq 1 22`;do
let k=${k}+1
if [ ${k} -eq ${SLURM_ARRAY_TASK_ID} ]; then


ldsc=/net/mulan/home/yasheng/comparisonProject/program/ldsc/ldsc.py
cd /net/mulan/home/yasheng/comparisonProject/program/ldsc
source activate ldsc

mkdir -p ${fbStr}reference/ukb/ldsc

${ldsc} --out ${fbStr}reference/ukb/ldsc/${chr} \
		--bfile ${fbStr}reference/ukb/geno/chr${chr} \
		--l2  --ld-wind-kb 1000.0 

fi
done

