#!/bin/bash


#SBATCH --partition=mulan,main
#SBATCH --time=12:00:00
#SBATCH --job-name=gemma
#SBATCH --mem=2G
#SBATCH --cpus-per-task=1
#SBATCH --array=1-1800%400
#SBATCH --output=/net/mulan/home/fredboe/research/ukb-intervals-sims/cluster_outputs/gemma_%a.out
#SBATCH --error=/net/mulan/home/fredboe/research/ukb-intervals-sims/cluster_outputs/gemma_%a.err


scenarios=( I II III IV)
distributions=( laplace normal scaledt)
hsqs=(0.1 0.2 0.5)
gemma=/net/fantasia/home/borang/software/gemma-0.98.3-linux-static


# 
let k=0
#let pc_ctr=0
type=ukb
let chr=1


for scenario in ${scenarios[@]}; do
  for distribution in ${distributions[@]}; do
    for hsq in ${hsqs[@]}; do


for p in `seq 1 10`; do
for fold in `seq 1 5`; do

let k=${k}+1
if [ ${k} -eq ${SLURM_ARRAY_TASK_ID} ]
then

dir=~/research/ukb-intervals-sims/dat-quant-136k/gemma/scenario${scenario}/${distribution}/hsq${hsq}
bfile=${dir}/chr${chr}
summ=summary_${type}_pheno${p}_scenario${scenario}_${distribution}_hsq${hsq}_fold${fold}_chr${chr}
let col=(${p}-1)*5+${fold}
echo ${col}


cd ${dir}
file=output/${summ}.assoc.txt
#if [ ! -f "$file" ]; then # check if file doesn't exist
${gemma} -bfile ${bfile} -notsnp -lm 1 -n ${col} -o ${summ}
#fi
fi

done 
done
done
done
done