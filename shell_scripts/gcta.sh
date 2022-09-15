#!/bin/bash


#SBATCH --partition=mulan,main
#SBATCH --time=12:00:00
#SBATCH --job-name=gcta
#SBATCH --mem=32G
#SBATCH --cpus-per-task=1
#SBATCH --array=1-36
#SBATCH --output=/net/mulan/home/fredboe/research/ukb-intervals-sims/cluster_outputs/gcta_%a.out
#SBATCH --error=/net/mulan/home/fredboe/research/ukb-intervals-sims/cluster_outputs/gcta_%a.err


declare -a files
for file in ~/research/ukb-intervals-sims/hapmap3/snp_effects/*.txt
do
    files=("${files[@]}" "$file")
done

#for filename in "${files[@]}"
let k=0

for filename in "${files[@]}"; do
  let k=${k}+1
  if [ ${k} -eq ${SLURM_ARRAY_TASK_ID} ]; then

    # parse hsq from file name
    # first, cut with period as field separator
    #filename=~/research/ukb-intervals-sims/hapmap3/snp_effects/scenarioIII_laplace_hsq0.5.txt
    out_suffix=$( basename $filename )
    field2=$( echo "${filename}" | cut -d'.' -f2 )
#    echo "$field2"
    hsq=$( echo "scale=2 ; $field2 / 10" | bc )
#    echo ${hsq}


    gcta64 --bfile ../hapmap3/plink_files_for_sims/chr1 \
            --simu-qt \
            --simu-causal-loci ${filename} \
            --simu-hsq ${hsq} \
            --simu-rep 10 \
            --out ../hapmap3/sim_traits/sims_${out_suffix}
    fi
  done
