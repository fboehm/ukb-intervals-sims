#!/bin/bash

#SBATCH --partition=mulan,main
#SBATCH --time=1-00:00:00
#SBATCH --job-name=SCT_u
#SBATCH --mem=60G
#SBATCH --cpus-per-task=5
#SBATCH --array=1-450%25
#SBATCH --output=/net/mulan/home/fredboe/research/ukb-intervals-sims/cluster_outputs/01_CT_ukb_b_thread5_%a.out
#SBATCH --error=/net/mulan/home/fredboe/research/ukb-intervals-sims/cluster_outputs/01_CT_ukb_b_thread5_%a.err

scenarios=( I )
distributions=( laplace normal scaledt)
hsqs=(0.1 0.2 0.5)

let k=0

#compstr=/net/mulan/disk2/yasheng/comparisonProject/
CT=../Rscript/01_CT_SCT.R
ref=ukb
dat=continuous
thread=5

for scenario in ${scenarios[@]}; do
  for distribution in ${distributions[@]}; do
    for hsq in ${hsqs[@]}; do
for p in `seq 1 10`; do
for fold in `seq 1 5`; do

let k=${k}+1

if [ ${k} -eq ${SLURM_ARRAY_TASK_ID} ]
then

if [[ "$dat" == "continuous" ]]
then
echo continuous
summ=../dat-quant/gemma/scenario${scenario}/${distribution}/hsq${hsq}/output/summary_ukb_pheno${p}_scenario${scenario}_${distribution}_hsq${hsq}_fold${fold}_chr1
#summ=${compstr}05_internal_c/pheno${p}/output/summary_ukb_cross${cross}
pathCT=../dat-quant/CT/scenario${scenario}/${distribution}/hsq${hsq}/
pathSCT=../dat-quant/SCT/scenario${scenario}/${distribution}/hsq${hsq}/
else
echo binary
summ=${compstr}06_internal_b/pheno${p}/output/summary_ukb_cross${cross}
path=${compstr}06_internal_b/pheno${p}/
fi

# cat ${summ}_chr*.assoc.txt > ${summ}.assoc.txt
# sed -i '/chr/d' ${summ}.assoc.txt
mkdir -p ${pathCT}
mkdir -p ${pathSCT}
file=${pathCT}esteff_ukb_pheno${p}_cross${fold}.txt

#esttime=${compstr}01_time_file/01_SCT_CT_ukb_${dat}_pheno${p}_cross${cross}_thread${thread}.tm
#time /usr/bin/time -v -o ${esttime} 

if [ ! -f "$file" ]; then # check if file doesn't exist

Rscript ${CT} --summ ${summ}.assoc.txt \
--pathCT ${pathCT} --pathSCT ${pathSCT} --pheno ${p} --cross ${fold} --reftype ${ref} --dat ${dat} --thread ${thread} \
--scenario ${scenario} --distribution ${distribution} --hsq ${hsq}

fi
# rm ${path}CT/summary_cross${cross}.bk
# rm ${path}CT/summary_cross${cross}.rds

fi
done
done
done
done
done

