#!/bin/bash

#SBATCH --partition=nomosix
#SBATCH --time=10:00:00
#SBATCH --job-name=sblup
#SBATCH --mem=63G
#SBATCH --cpus-per-task=5

#SBATCH --array=1-110%50
#SBATCH --output=/net/mulan/disk2/yasheng/comparisonProject/00_cluster_file/05_sblup_ukb_thread5_%a.out
#SBATCH --error=/net/mulan/disk2/yasheng/comparisonProject/00_cluster_file/05_sblup_ukb_thread5_%a.err

bash
let k=0
let ldr=200
let thread=5

gcta=/net/mulan/home/yasheng/comparisonProject/program/gcta_1.93.1beta/gcta64
compstr=/net/mulan/disk2/yasheng/comparisonProject/
dat=continuous

for p in 1 ;do
for cross in 1 2 3 4 5; do
for chr in `seq 1 22`; do

let k=${k}+1
if [ ${k} -eq ${SLURM_ARRAY_TASK_ID} ]; then

if [[ "$dat" == "continuous" ]]
then
echo continuous
herit=${compstr}05_internal_c/pheno${p}/herit/h2_ukb_cross${cross}.log
summ=${compstr}05_internal_c/pheno${p}/output/summary_ukb_cross${cross}_chr${chr}
est=/net/mulan/disk2/yasheng/comparisonProject-archive/05_internal_c/pheno${p}/sblup/esteff_ukb_cross${cross}_chr${chr}
fi

if [[ "$dat" == "binary" ]]
then
echo binary
herit=${compstr}06_internal_b/pheno${p}/herit/h2_ukb_cross${cross}.log
summ=${compstr}06_internal_b/pheno${p}/output/summary_ukb_cross${cross}_chr${chr}
est=${compstr}06_internal_b/pheno${p}/sblup/esteff_ukb_cross${cross}_chr${chr}
fi

## heritability
hstr=`sed -n '26p' ${herit}`
hse=`echo ${hstr#*:}`
h2=`echo ${hse%(*}`

## snp number
mstr=`sed -n '24p' ${herit}`
mstrr=`echo ${mstr#*,}`
m=`echo ${mstrr% S*}`
cojo=$(echo "${m}*(1/${h2}-1)" | bc -l)

## sblup estimation
awk '{print $2,$6,$7,$8,$9,$10,$11,$5}' ${summ}.assoc.txt > ${summ}.ma
sed -i '1i\SNP A1 A2 freq b se p N' ${summ}.ma
ref=${compstr}04_reference/ukb/geno/chr${chr}
esttime=${compstr}01_time_file/05_sblup_ukb_${dat}_pheno${p}_cross${cross}_chr${chr}_thread${thread}.tm
time /usr/bin/time -v -o ${esttime} ${gcta} --bfile ${ref} --chr ${chr} --cojo-file ${summ}.ma --cojo-sblup ${cojo} --cojo-wind ${ldr} --thread-num ${thread} --out ${est} 

## remove file
rm ${est}.*badsnps
rm ${summ}.ma
rm ${est}.log

fi
done
done
done
