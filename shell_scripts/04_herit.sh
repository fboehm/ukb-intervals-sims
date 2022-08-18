#!/bin/bash

#SBATCH --partition=mulan,main
#SBATCH --time=04:00:00
#SBATCH --job-name=herit
#SBATCH --mem-per-cpu=10G
#SBATCH --array=1-25
#SBATCH --output=/net/mulan/home/fredboe/research/ukb-intervals-sims/cluster_outputs/04_herit_%j_%a.out
#SBATCH --error=/net/mulan/home/fredboe/research/ukb-intervals-sims/cluster_outputs/04_herit_%j_%a.err

pc=0.1
hsq=0.1

let k=0

ldsc=/net/mulan/home/yasheng/comparisonProject/program/ldsc/ldsc.py
mkldsc=/net/mulan/home/yasheng/comparisonProject/code/02_method/04_mk_ldsc_summ.R
compstr=/net/mulan/disk2/yasheng/comparisonProject/
dat=c
reftype=ukb


for p in `seq 1 5`; do
for cross in `seq 1 5`; do
let k=${k}+1
if [ ${k} -eq ${SLURM_ARRAY_TASK_ID} ]; then

ref=/net/mulan/disk2/yasheng/comparisonProject/04_reference/${reftype}/ldsc/

if [[ "$dat" == "c" ]]
then
summ1=~/research/ukb-intervals-sims/dat/hsq${hsq}_pcausal${pc}/gemma/output/summary_ukb_pheno${p}_fold${cross}_chr
summ2=~/research/ukb-intervals-sims/dat/hsq${hsq}_pcausal${pc}/gemma/output/summary_ukb_pheno${p}_fold${cross}

#summ=${compstr}05_internal_c/pheno${p}/output/summary_${reftype}_cross${cross}
#h2=${compstr}05_internal_c/pheno${p}/herit/h2_${reftype}_cross${cross}
h2path=~/research/ukb-intervals-sims/dat/hsq${hsq}_pcausal${pc}/ldsc
h2=${h2path}/h2_${reftype}_fold${cross}_pheno${p}_hsq${hsq}_pcausal${pc}
else
summ=${compstr}06_internal_b/pheno${p}/output/summary_hm3_cross${cross}
h2=${compstr}06_internal_b/pheno${p}/herit/h2_${reftype}_cross${cross}
fi

## summary data for ldsc
cat ${summ1}*.assoc.txt > ${summ2}.assoc.txt
sed -i '/chr/d' ${summ2}.assoc.txt
Rscript ${mkldsc} --summgemma ${summ2}.assoc.txt --summldsc ${summ2}.ldsc

## heritability
#source activate /net/mulan/home/yasheng/py3/envs/ldsc
mkdir -p ${h2path}
python2 ${ldsc} --h2 ${summ2}.ldsc.gz --ref-ld-chr ${ref} --w-ld-chr ${ref} --out ${h2}
fi
done
done
done
