#!/bin/bash



#SBATCH --partition=mulan,main
#SBATCH --time=1-00:00:00
#SBATCH --job-name=subset-val-ver
#SBATCH --mem-per-cpu=2G
#SBATCH --array=1-22
#SBATCH --output=/net/mulan/home/fredboe/research/ukb-intervals-sims/cluster_outputs/process_dat-05_get_subsample_validation_%j_%a.out
#SBATCH --error=/net/mulan/home/fredboe/research/ukb-intervals-sims/cluster_outputs/process_dat-05_get_subsample_validation_%j_%a.err




compStr=/net/mulan/disk2/yasheng/comparisonProject/
fbStr=~/research/ukb-intervals-sims/dat/

# NOTE THAT THE *SIMULATED* DATA ALL HAVE THE SUBJECTS SUBSETTED IN TEH SAME WAY
# IE, THEY HAVE THE SAME VALIDATION, VERIFICATION, AND REFERENCE SUBSETS
# ukb subsample
let k=0
for chr in `seq 1 22`;do
let k=${k}+1
if [ ${k} -eq ${SLURM_ARRAY_TASK_ID} ]; then

bfileAll=/net/mulan/disk2/yasheng/predictionProject/plink_file/ukb/chr${chr}

## validation sample
bfileSub=${fbStr}validation/geno/allchr${chr}
bfileSubP=${fbStr}validation/geno/chr${chr}
mkdir -p ${fbStr}validation/geno
idxSub=${fbStr}validation-ids.txt
plink-1.9 --silent --bfile ${bfileAll} --keep ${idxSub} --make-bed --out ${bfileSub}
plink-1.9 --silent --bfile ${bfileSub} --maf 0.01 --make-bed --out ${bfileSubP}
rm ${bfileSub}.log
rm ${bfileSub}.bed
rm ${bfileSub}.bim
rm ${bfileSub}.fam
rm ${bfileSubP}.log

# impute for validation
geno_impute=${compStr}code/01_process_dat/05_geno_imputation.R
input=${bfileSubP}
output=${fbStr}validation/geno/impute/chr${chr}
mkdir -p ${fbStr}validation/geno/impute
Rscript ${geno_impute} --plinkin ${input} --plinkout ${output}

### verification data

bfileSub=${fbStr}verification/geno/allchr${chr}
bfileSubP=${fbStr}verification/geno/chr${chr}
mkdir -p ${fbStr}verification/geno
idxSub=${fbStr}verification-ids.txt
plink-1.9 --silent --bfile ${bfileAll} --keep ${idxSub} --make-bed --out ${bfileSub}
plink-1.9 --silent --bfile ${bfileSub} --maf 0.01 --make-bed --out ${bfileSubP}
rm ${bfileSub}.log
rm ${bfileSub}.bed
rm ${bfileSub}.bim
rm ${bfileSub}.fam
rm ${bfileSubP}.log

# impute for verification
geno_impute=${compStr}code/01_process_dat/05_geno_imputation.R
input=${bfileSubP}
output=${fbStr}verification/geno/impute/chr${chr}
mkdir -p ${fbStr}verification/geno/impute

Rscript ${geno_impute} --plinkin ${input} --plinkout ${output}

fi
done


