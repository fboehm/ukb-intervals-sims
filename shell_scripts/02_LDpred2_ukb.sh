#!/bin/bash

#SBATCH --partition=mulan,main
#SBATCH --time=1-00:00:00
#SBATCH --job-name=LDpred2_ukb
#SBATCH --mem=30G
#SBATCH --cpus-per-task=5

#SBATCH --array=1-550%110
#SBATCH --output=/net/mulan/home/fredboe/research/ukb-intervals-sims/cluster_outputs/02_LDpred2_ukb_thread5_%j_%a.out
#SBATCH --error=/net/mulan/home/fredboe/research/ukb-intervals-sims/cluster_outputs/02_LDpred2_ukb_thread5_%j_%a.err


hsq=0.1
pc=0.1

let k=0
let thread=${SLURM_CPUS_PER_TASK}

compstr=/net/mulan/disk2/yasheng/comparisonProject/
LDpred2=../Rscript/02_LDpred2.R
ref=ukb
dat=continuous

for model in LDpred2-auto;do
for p in `seq 1 5`; do
for cross in 1 2 3 4 5; do
for chr in `seq 1 22`;do
let k=${k}+1
if [ ${k} -eq ${SLURM_ARRAY_TASK_ID} ]; then

if [[ "$dat" == "continuous" ]]
then
summ=~/research/ukb-intervals-sims/dat/hsq${hsq}_pcausal${pc}/gemma/output/summary_ukb_fold${cross}_chr
LDpred2Path=~/research/ukb-intervals-sims/dat/hsq${hsq}_pcausal${pc}/LDpred2/
else
summ=/net/mulan/disk2/yasheng/comparisonProject/06_internal_b/pheno${p}/output/summary_ukb_pheno${p}_fold${cross}_chr
LDpred2Path=/net/mulan/disk2/yasheng/comparisonProject/06_internal_b/pheno${p}/LDpred2/
fi
mkdir -p ${LDpred2Path}
esttime=~/research/ukb-intervals-sims/cluster_outputs/02_${model}_ukb_${dat}_pheno${p}_cross${cross}_chr${chr}_thread${thread}.tm
time /usr/bin/time -v -o ${esttime} Rscript ${LDpred2} --summ ${summ} --LDpred2Path ${LDpred2Path} --pheno ${p} --cross ${cross} --dat ${dat} --reftype ${ref} --thread ${thread} --chr ${chr} --model ${model}

fi
done
done
done
done
