#!/bin/bash

#SBATCH --partition=mulan,main
#SBATCH --time=3-00:00:00
#SBATCH --job-name=cvplus
#SBATCH --mem=16G
#SBATCH --output=/net/mulan/home/fredboe/research/ukb-intervals-sims/cluster_outputs/render-rmd_%j_%a.out
#SBATCH --error=/net/mulan/home/fredboe/research/ukb-intervals-sims/cluster_outputs/render-rmd_%j_%a.err

Rscript -e 'rmarkdown::render(input = "~/research/ukb-intervals-sims/Rmd/cv-plus-simulations.Rmd",
                              output_file = "~/research/ukb-intervals-sims/results/DBSLMM/cv-plus-simulations-alpha0.1-hsq0.1-pcausal0.1.html",
                              params = list(hsq=0.1, pcausal=0.1))'
