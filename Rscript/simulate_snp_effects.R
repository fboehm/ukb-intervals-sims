#!/usr/bin/env Rscript

library(magrittr)

#bim_file <- "~/research/ukb-intervals/dat/plink_files/ukb/chr22.bim"
hsq <- snakemake@config[["hsq"]]
p_causal <- snakemake@config[["pcausal"]]
bim_file <- snakemake@input[[1]]
bim_tib <- vroom::vroom(bim_file, col_names = FALSE) %>%
  dplyr::mutate(snp_index = 1:nrow(.)) %>%
  dplyr::rename(chromosome = X1, snp_id = X2, pos_cm = X3, pos_bp = X4, a1 = X5, a2 = X6)

# set params
m <- nrow(bim_tib)
p_causal <- args[2]
hsq <- args[1]
set.seed(2022-05-23)

bim_tib2 <- bim_tib %>%
  dplyr::mutate(causal = rbinom(n = m, size = 1, prob = p_causal)) %>%
  dplyr::mutate(effect = causal * rnorm(n = m, mean = 0, sd = sqrt(hsq / (m * p_causal))))
# make causal.snplist for GCTA
bim_tib2 %>%
  dplyr::select(snp_id, effect) %>%
  vroom::vroom_write(file = snakemake@output[[1]], col_names = FALSE)
  #vroom::vroom_write(file = "../dat/simulations-ding/snp_effects_Chr22_hsq0.2_pcausal0.1.txt",
  #                   col_names = FALSE)
