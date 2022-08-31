#!/bin/bash

#SBATCH --partition=mulan,main
#SBATCH --time=1-00:00:00
#SBATCH --job-name=lassosum_ukb
#SBATCH --mem=40G
#SBATCH --cpus-per-task=1

#SBATCH --array=1-125%25
#SBATCH --output=/net/mulan/home/fredboe/research/ukb-intervals-sims/cluster_outputs/03_lassosum_ukb_thread1_%j_%a.out
#SBATCH --error=/net/mulan/home/fredboe/research/ukb-intervals-sims/cluster_outputs/03_lassosum_ukb_thread1_%j_%a.err


hsq=0.1
pc=0.1

let k=0

ref=ukb
dat=continuous
#thread=4
thread=1
compstr=/net/mulan/disk2/yasheng/comparisonProject/
lassosum=../Rscript/03_lassosum.R


for p in `seq 1 5`; do
for cross in 1 2 3 4 5; do

let k=${k}+1
if [ ${k} -eq ${SLURM_ARRAY_TASK_ID} ]
then

if [[ "$dat" == "continuous" ]]
then
echo continuous 
summ=~/research/ukb-intervals-sims/dat/hsq${hsq}_pcausal${pc}/gemma/output/summary_ukb_pheno${p}_fold${cross}
path=~/research/ukb-intervals-sims/dat/hsq${hsq}_pcausal${pc}/
else
echo binary
summ=${compstr}06_internal_b/pheno${p}/output/summary_ukb_cross${cross}
path=${compstr}06_internal_b/pheno${p}/lassosum/
fi

cat ${summ}_chr*.assoc.txt > ${summ}.assoc.txt
sed -i '/chr/d' ${summ}.assoc.txt

nobs=`sed -n "2p" ${summ}.assoc.txt | awk '{print $5}'`
nmis=`sed -n "2p" ${summ}.assoc.txt | awk '{print $4}'`
n=$(echo "${nobs}+${nmis}" | bc -l)

mkdir -p ${path}
esttime=~/research/ukb-intervals-sims/cluster_outputs/03_lassosum_ukb_${dat}_pheno${p}_cross${cross}_thread${thread}.tm
time /usr/bin/time -v -o ${esttime} Rscript ${lassosum} --summ ${summ}.assoc.txt --dat ${dat} --pheno ${p} --n ${n} --thread ${thread} --cross ${cross} --reftype ${ref} --path ${path}
fi

done
done
