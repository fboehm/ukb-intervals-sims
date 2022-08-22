#!/bin/bash


#SBATCH --partition=mulan,main
#SBATCH --time=18:00:00
#SBATCH --job-name=gcta
#SBATCH --mem=32G
#SBATCH --cpus-per-task=1
#SBATCH --output=/net/mulan/home/fredboe/research/ukb-intervals-sims/cluster_outputs/gcta_%j_%a.out
#SBATCH --error=/net/mulan/home/fredboe/research/ukb-intervals-sims/cluster_outputs/gcta_%j_%a.err

hsq=0.5
pc=0.1

mkdir -p ../dat/hsq${hsq}_pcausal${pc}/sim_traits
    
gcta64 --bfile ../dat/chr22 \
        --simu-qt \
        --simu-causal-loci ../dat/hsq${hsq}_pcausal${pc}/snp_effects_Chr22_hsq${hsq}_pcausal${pc}.txt \
        --simu-hsq ${hsq} \
        --simu-rep 5 \
        --out ../dat/hsq${hsq}_pcausal${pc}/sim_traits/sims_Chr22_hsq${hsq}_pcausal${pc}
            
