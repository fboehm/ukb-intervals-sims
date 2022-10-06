#!/bin/bash

#SBATCH --partition=mulan,main
#SBATCH --time=10:00
#SBATCH --job-name=calc_pgs
#SBATCH --cpus-per-task=1
#SBATCH --mem=32G
#SBATCH --array=1-360%50
#SBATCH --output=/net/mulan/home/fredboe/research/ukb-intervals-sims/cluster_outputs/calc_pgs_%a.out
#SBATCH --error=/net/mulan/home/fredboe/research/ukb-intervals-sims/cluster_outputs/calc_pgs_%a.err

scenarios=( I II III IV)
distributions=( laplace normal scaledt)
hsqs=(0.1 0.2 0.5)
let k=0

for scenario in ${scenarios[@]}; do
  for distribution in ${distributions[@]}; do
    for hsq in ${hsqs[@]}; do
      for p in `seq 1 10`; do
        let k=${k}+1
        if [ ${k} -eq ${SLURM_ARRAY_TASK_ID} ]; then
          Rscript --vanilla ~/research/ukb-intervals-sims/Rscript/calc_pgs.R ${scenario} ${distribution} ${hsq} ${p}      
        fi
      done
    done
  done
done