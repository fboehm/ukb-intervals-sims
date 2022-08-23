#!/bin/bash


#SBATCH --partition=mulan,main
#SBATCH --time=3-00:00:00
#SBATCH --job-name=gemma
#SBATCH --mem=2G
#SBATCH --cpus-per-task=1
#SBATCH --array=1-550%70
#SBATCH --output=/net/mulan/home/fredboe/research/ukb-intervals-sims/cluster_outputs/gemma_%j_%a.out
#SBATCH --error=/net/mulan/home/fredboe/research/ukb-intervals-sims/cluster_outputs/gemma_%j_%a.err


hsq=0.5
pc=0.1


# 
let k=0
#let pc_ctr=0
type=ukb
gemma=/net/mulan/home/yasheng/comparisonProject/program/gemma-0.98.1-linux-static

for p in `seq 1 5`
do
for fold in 1 2 3 4 5
do 
for chr in `seq 1 22`
do

let k=${k}+1
if [ ${k} -eq ${SLURM_ARRAY_TASK_ID} ]
then

bfile=~/research/ukb-intervals-sims/dat/hsq${hsq}_pcausal${pc}/gemma/chr${chr}
summ=summary_${type}_pheno${p}_fold${fold}_chr${chr}

let col=(${p}-1)\*5+${fold}
echo ${col}
outpath=~/research/ukb-intervals-sims/dat/hsq${hsq}_pcausal${pc}/gemma
cd ${outpath}
file=output/${summ}.assoc.txt
if [ ! -f "$file" ]; then # check if file doesn't exist
${gemma} -bfile ${bfile} -notsnp -lm 1 -n ${col} -o ${summ}
fi
fi

done 
done
done
