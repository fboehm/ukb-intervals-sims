#!/bin/bash

#SBATCH --partition=mulan,main
#SBATCH --time=3-00:00:00
#SBATCH --job-name=DBSLMM
#SBATCH --mem=12G
#SBATCH --cpus-per-task=5

#SBATCH --array=1-360%30
#SBATCH --output=/net/mulan/home/fredboe/research/ukb-intervals-sims/cluster_outputs/06_DBSLMM-tuning_sims_c_%j_%a.out
#SBATCH --error=/net/mulan/home/fredboe/research/ukb-intervals-sims/cluster_outputs/06_DBSLMM-tuning_sims_c_%j_%a.err

scenarios=( I II III IV)
distributions=( laplace normal scaledt)
hsqs=(0.1 0.2 0.5)


# 
let k=0
#let pc_ctr=0
chr=1


for scenario in ${scenarios[@]}; do
  for distribution in ${distributions[@]}; do
    for hsq in ${hsqs[@]}; do
      for p in `seq 1 10`; do

let k=${k}+1
if [ ${k} -eq ${SLURM_ARRAY_TASK_ID} ]
then

let thread=5

dat=continuous
type=t

fbStr=~/research/ukb-intervals-sims/

compstr=/net/mulan/disk2/yasheng/comparisonProject/
plink=/usr/cluster/bin/plink-1.9
DBSLMM=06_DBSLMM_script.sh
DBSLMMpath=/net/mulan/home/yasheng/predictionProject/code/ # path to dbslmm C++ compiled
blockf=${compstr}LDblock_EUR/chr
#ref=${compstr}04_reference/ukb/geno/chr
ref=${fbStr}dat-quant/reference/chr

let k=${k}+1
if [ ${k} -eq ${SLURM_ARRAY_TASK_ID} ]; then

val=~/research/ukb-intervals-sims/dat-quant/validation/chr

if [[ "$dat" == "continuous" ]]
then
phenoVal=~/research/ukb-intervals-sims/dat-quant/validation/scenario${scenario}_distribution${distribution}_hsq${hsq}_replicate${p}.txt
index=r2
else
phenoVal=${compstr}03_subsample/${dat}/pheno${p}/val/ukb/02_pheno_b.txt
index=auc
fi

## input
if [[ "$dat" == "continuous" ]]
then
herit=${hsq}
#summ=${compstr}05_internal_c/pheno${p}/output/summary_ukb_cross${cross}_chr
summ=~/research/ukb-intervals-sims/dat/hsq${hsq}_pcausal${pcausal}/gemma/output/summary_ukb_pheno${p}_fold${cross}_chr
#outPath=/net/mulan/disk2/yasheng/comparisonProject/05_internal_c/pheno${p}/DBSLMM/
outPath=~/research/ukb-intervals-sims/dat-quant/DBSLMM/scenario${scenario}/${distribution}/hsq${hsq}
else
herit=${compstr}06_internal_b/pheno${p}/herit/h2_ukb_cross${cross}.log
summ=${compstr}06_internal_b/pheno${p}/output/summary_ukb_cross${cross}_chr
outPath=/net/mulan/disk2/yasheng/comparisonProject/06_internal_b/pheno${p}/DBSLMM/
cov=${compstr}03_subsample/${dat}/pheno${p}/val/ukb/03_cov_eff.txt
fi

mkdir -p ${outPath}
## DBSLMM
esttime=~/research/ukb-intervals-sims/cluster_outputs/06_DBSLMM_ukb_scenario${scenario}_distribution${distribution}_hsq${hsq}_pheno${p}.tm
if [[ "$dat" == "continuous" ]]
then
time /usr/bin/time -v -o ${esttime} sh ${DBSLMM} -D ${DBSLMMpath} -p ${plink} -B ${blockf} -s ${summ} -m DBSLMM\
             -H ${herit} -G ${val} -R ${ref} -P ${phenoVal}\
             -l 1 -T ${type} -i ${index} -t ${thread} -o ${outPath} #-C ${chr} -f ${h2f} -h ${pth}
else 
# time /usr/bin/time -v -o ${esttime} 
sh ${DBSLMM} -D ${DBSLMMpath} -p ${plink} -B ${blockf} -s ${summ}  -m DBSLMM\
             -H ${herit} -G ${val} -R ${ref} -P ${phenoVal}\
             -l 1 -T ${type} -c ${cov} -i ${index} -t ${thread} -o ${outPath} #-C ${chr} -f ${h2f} -h ${pth}
fi


# for chr in `seq 1 22`
# do 
# gzip ${summ}${chr}.assoc.txt
# done

fi
done
done
done
done


