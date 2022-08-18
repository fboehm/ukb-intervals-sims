#!/usr/bin/env Rscript

hsq <- 0.1
pcausal <- 0.1

library(magrittr)


bim_file <- "../dat/chr22.bim"
bim_tib <- vroom::vroom(bim_file, col_names = FALSE) %>%
  dplyr::mutate(snp_index = 1:nrow(.)) %>%
  dplyr::rename(chromosome = X1, snp_id = X2, pos_cm = X3, pos_bp = X4, a1 = X5, a2 = X6)

# set params
m <- nrow(bim_tib)
set.seed(2022-05-23)

bim_tib2 <- bim_tib %>%
  dplyr::mutate(causal = rbinom(n = m, size = 1, prob = pcausal)) %>%
  dplyr::mutate(effect = causal * rnorm(n = m, mean = 0, sd = sqrt(hsq / (m * pcausal))))
# make causal.snplist for GCTA
bim_tib2 %>%
  dplyr::select(snp_id, effect) %>%
  vroom::vroom_write(file = paste0("../dat/hsq", hsq, "_pcausal", pcausal, "/snp_effects_Chr22_hsq", hsq, "_pcausal", pcausal, ".txt"), col_names = FALSE)
