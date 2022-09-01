#!/bin/bash

#SBATCH --partition=mulan,main
#SBATCH --time=3-00:00:00
#SBATCH --job-name=DBSLMM
#SBATCH --mem=12G
#SBATCH --cpus-per-task=5

#SBATCH --array=1-25%5
#SBATCH --output=/net/mulan/home/fredboe/research/ukb-intervals-sims/cluster_outputs/06_DBSLMM-tuning_sims_c_%j_%a.out
#SBATCH --error=/net/mulan/home/fredboe/research/ukb-intervals-sims/cluster_outputs/06_DBSLMM-tuning_sims_c_%j_%a.err

hsq=0.1
pcausal=0.1


let k=0
let thread=5

dat=continuous
type=t

compstr=/net/mulan/disk2/yasheng/comparisonProject/
plink=/usr/cluster/bin/plink-1.9
DBSLMM=${compstr}code/02_method/06_DBSLMM_script.sh
#DBSLMM=06_DBSLMM_script.sh
DBSLMMpath=/net/mulan/home/yasheng/predictionProject/code/
blockf=${compstr}LDblock_EUR/chr
#ref=${compstr}04_reference/ukb/geno/chr


for p in `seq 1 5`; do
for cross in 1 2 3 4 5; do
let k=${k}+1
if [ ${k} -eq ${SLURM_ARRAY_TASK_ID} ]; then

# phenoMiss=/net/mulan/disk2/yasheng/comparisonProject/code/02_method/DBSLMM_miss/pheno.txt
# crossMiss=/net/mulan/disk2/yasheng/comparisonProject/code/02_method/DBSLMM_miss/cross.txt
# chrMiss=/net/mulan/disk2/yasheng/comparisonProject/code/02_method/DBSLMM_miss/chr.txt
# h2fMiss=/net/mulan/disk2/yasheng/comparisonProject/code/02_method/DBSLMM_miss/h2f.txt
# pthMiss=/net/mulan/disk2/yasheng/comparisonProject/code/02_method/DBSLMM_miss/pth.txt
# for iter in `seq 1 5`
# do
# let k=${k}+1
# if [ ${k} -eq ${SLURM_ARRAY_TASK_ID} ]
# then
# p=`head -n ${iter} ${phenoMiss} | tail -n 1`
# cross=`head -n ${iter} ${crossMiss} | tail -n 1`
# chr=`head -n ${iter} ${chrMiss} | tail -n 1`
# h2f=`head -n ${iter} ${h2fMiss} | tail -n 1`
# pth=`head -n ${iter} ${pthMiss} | tail -n 1`
# echo pheno${p}_cross${cross}_chr${chr}_h2f${h2f}_pth${pth}

# 
#val=${compstr}03_subsample/${dat}/pheno${p}/val/ukb/impute_inter/chr
val=~/research/ukb-intervals-sims/dat/hsq${hsq}_pcausal${pcausal}/validation/chr

if [[ "$dat" == "continuous" ]]
then
# phenoVal=${compstr}03_subsample/${dat}/pheno${p}/02_pheno_c.txt
#phenoVal=${compstr}/03_subsample/${dat}/pheno${p}/val/ukb/02_pheno_c.txt
phenoVal=~/research/ukb-intervals-sims/dat/hsq${hsq}_pcausal${pcausal}/validation/pheno${p}_hsq${hsq}_pcausal${pcausal}.txt
index=r2
else
phenoVal=${compstr}03_subsample/${dat}/pheno${p}/val/ukb/02_pheno_b.txt
index=auc
fi

## input
if [[ "$dat" == "continuous" ]]
then
#herit=${compstr}05_internal_c/pheno${p}/herit/h2_ukb_cross${cross}.log
herit=~/research/ukb-intervals-sims/dat/hsq${hsq}_pcausal${pcausal}/ldsc/h2_ukb_fold${cross}_pheno${p}_hsq${hsq}_pcausal${pcausal}.log
#summ=${compstr}05_internal_c/pheno${p}/output/summary_ukb_cross${cross}_chr
summ=~/research/ukb-intervals-sims/dat/hsq${hsq}_pcausal${pcausal}/gemma/output/summary_ukb_pheno${p}_fold${cross}_chr
#outPath=/net/mulan/disk2/yasheng/comparisonProject/05_internal_c/pheno${p}/DBSLMM/
outPath=~/research/ukb-intervals-sims/dat/hsq${hsq}_pcausal${pcausal}/DBSLMM/
else
herit=${compstr}06_internal_b/pheno${p}/herit/h2_ukb_cross${cross}.log
summ=${compstr}06_internal_b/pheno${p}/output/summary_ukb_cross${cross}_chr
outPath=/net/mulan/disk2/yasheng/comparisonProject/06_internal_b/pheno${p}/DBSLMM/
cov=${compstr}03_subsample/${dat}/pheno${p}/val/ukb/03_cov_eff.txt
fi

mkdir -p ${outPath}
## DBSLMM
esttime=~/research/ukb-intervals-sims/cluster_outputs/06_DBSLMM_ukb_${dat}_pheno${p}_cross${cross}_hsq${hsq}_pcausal${pcausal}_thread${thread}.tm
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



