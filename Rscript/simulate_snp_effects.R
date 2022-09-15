#!/usr/bin/env Rscript


library(magrittr)


bim_file <- "../hapmap3/hm3_ukb/chr1.bim"
bim_tib <- vroom::vroom(bim_file, col_names = FALSE) %>%
  dplyr::mutate(snp_index = 1:nrow(.)) %>%
  dplyr::rename(chromosome = X1, snp_id = X2, pos_cm = X3, pos_bp = X4, a1 = X5, a2 = X6)

# set params
m <- nrow(bim_tib)
set.seed(2022-05-23)

## SCENARIO I
hsq_vec <- c(0.1, 0.2, 0.5)
# normal distribution
for (hsq in hsq_vec){
  outfile <- paste0("../hapmap3/snp_effects/scenarioI_normal_hsq", hsq, ".txt")
  bim_tib %>%
    dplyr::mutate(effect = rnorm(n = m, mean = 0, sd = sqrt(hsq / m))) %>%
    dplyr::select(snp_id, effect) %>%
    vroom::vroom_write(file = outfile, col_names = FALSE)
}
# scaled t distribution with 4 df
# The scaled, shifted t distribution has mean mean and variance sd^2 * df/(df-2)
# here, df / (df -2) has value 2
# To get the sum of m independent rvs variances to be hsq, we need an average of hsq / m 
# ie, 2sd^2 = hsq / m; so sd = sqrt(hsq / (2m))
for (hsq in hsq_vec){
  outfile <- paste0("../hapmap3/snp_effects/scenarioI_scaledt_hsq", hsq, ".txt")
  bim_tib %>%
    dplyr::mutate(effect = metRology::rt.scaled(n = m, df = 4, mean = 0, sd = sqrt(hsq / (2 * m)))) %>%
    dplyr::select(snp_id, effect) %>%
    vroom::vroom_write(file = outfile, col_names = FALSE)
}
# laplace distribution
# var = 2 * scale^2
for (hsq in hsq_vec){
  outfile <- paste0("../hapmap3/snp_effects/scenarioI_laplace_hsq", hsq, ".txt")
  bim_tib %>%
    dplyr::mutate(effect = VGAM::rlaplace(n = m, location = 0, scale = sqrt(hsq / (2 * m)))) %>%
    dplyr::select(snp_id, effect) %>%
    vroom::vroom_write(file = outfile, col_names = FALSE)
}

## SCENARIO II
pcausal <- 0.001
# normal distribution
for (hsq in hsq_vec){
  outfile <- paste0("../hapmap3/snp_effects/scenarioII_normal_hsq", hsq, ".txt")
  bim_tib %>%
    dplyr::mutate(causal = rbinom(n = m, size = 1, prob = pcausal)) %>%
    dplyr::mutate(effect = causal * rnorm(n = m, mean = 0, sd = sqrt(hsq / sum(causal)))) %>%
    dplyr::select(snp_id, effect) %>%
    vroom::vroom_write(file = outfile, col_names = FALSE)
}
# scaled t distribution with 4 df
# The scaled, shifted t distribution has mean mean and variance sd^2 * df/(df-2)
# here, df / (df -2) has value 2
# To get the sum of m independent rvs variances to be hsq, we need an average of hsq / m 
# ie, 2sd^2 = hsq / m; so sd = sqrt(hsq / (2m))
for (hsq in hsq_vec){
  outfile <- paste0("../hapmap3/snp_effects/scenarioII_scaledt_hsq", hsq, ".txt")
  bim_tib %>%
    dplyr::mutate(causal = rbinom(n = m, size = 1, prob = pcausal)) %>%
    dplyr::mutate(effect = causal * metRology::rt.scaled(n = m, df = 4, mean = 0, sd = sqrt(hsq / (2 * sum(causal))))) %>%
    dplyr::select(snp_id, effect) %>%
    vroom::vroom_write(file = outfile, col_names = FALSE)
}
# laplace distribution
# var = 2 * scale^2
for (hsq in hsq_vec){
  outfile <- paste0("../hapmap3/snp_effects/scenarioII_laplace_hsq", hsq, ".txt")
  bim_tib %>%
    dplyr::mutate(causal = rbinom(n = m, size = 1, prob = pcausal)) %>%
    dplyr::mutate(effect = causal * VGAM::rlaplace(n = m, location = 0, scale = sqrt(hsq / (2 * sum(causal))))) %>%
    dplyr::select(snp_id, effect) %>%
    vroom::vroom_write(file = outfile, col_names = FALSE)
}

## Scenario III
pge <- 0.2 # proportion of genetic variance due to large effect SNPs
plarge <- 0.001 # proportion of snps with large effects
# normal distribution
for (hsq in hsq_vec){
  outfile <- paste0("../hapmap3/snp_effects/scenarioIII_normal_hsq", hsq, ".txt")
  bim_tib %>%
    dplyr::mutate(large = rbinom(n = m, size = 1, prob = plarge)) %>%
    dplyr::mutate(effect_large = rnorm(n = m, mean = 0, sd = sqrt(pge * hsq / sum(large)))) %>%
    dplyr::mutate(effect_small = rnorm(n = m, mean = 0, sd = sqrt((1 - pge) * hsq / (m - sum(large))))) %>%
    dplyr::mutate(effect = large * effect_large + (1 - large) * effect_small) %>%
    dplyr::select(snp_id, effect) %>%
    vroom::vroom_write(file = outfile, col_names = FALSE)
}
# scaled t distribution with 4 df
# The scaled, shifted t distribution has mean mean and variance sd^2 * df/(df-2)
# here, df / (df -2) has value 2
# To get the sum of m independent rvs variances to be hsq, we need an average of hsq / m 
# ie, 2sd^2 = hsq / m; so sd = sqrt(hsq / (2m))
for (hsq in hsq_vec){
  outfile <- paste0("../hapmap3/snp_effects/scenarioIII_scaledt_hsq", hsq, ".txt")
  bim_tib %>%
    dplyr::mutate(large = rbinom(n = m, size = 1, prob = plarge)) %>%
    dplyr::mutate(effect_large = metRology::rt.scaled(n = m, df = 4, mean = 0, sd = sqrt(pge * hsq / (2 * sum(large))))) %>%
    dplyr::mutate(effect_small = metRology::rt.scaled(n = m, df = 4, mean = 0, sd = sqrt((1 - pge) * hsq / (2 * (m - sum(large)))))) %>%
    dplyr::mutate(effect = large * effect_large + (1 - large) * effect_small) %>%
    dplyr::select(snp_id, effect) %>%
    vroom::vroom_write(file = outfile, col_names = FALSE)
}
# laplace distribution
# var = 2 * scale^2
for (hsq in hsq_vec){
  outfile <- paste0("../hapmap3/snp_effects/scenarioIII_laplace_hsq", hsq, ".txt")
  bim_tib %>%
    dplyr::mutate(large = rbinom(n = m, size = 1, prob = plarge)) %>%
    dplyr::mutate(effect_large = VGAM::rlaplace(n = m, location = 0, scale = sqrt(pge * hsq / (2 * sum(large))))) %>%
    dplyr::mutate(effect_small = VGAM::rlaplace(n = m, location = 0, scale = sqrt((1 - pge) * hsq / (2 * (m - sum(large)))))) %>%
    dplyr::mutate(effect = large * effect_large + (1 - large) * effect_small) %>%
    dplyr::select(snp_id, effect) %>%
    vroom::vroom_write(file = outfile, col_names = FALSE)
}

## Scenario IV
pge <- 0.5 # proportion of genetic variance due to large effect SNPs
plarge <- 0.001 # proportion of snps with large effects
# normal distribution
for (hsq in hsq_vec){
  outfile <- paste0("../hapmap3/snp_effects/scenarioIV_normal_hsq", hsq, ".txt")
  bim_tib %>%
    dplyr::mutate(large = rbinom(n = m, size = 1, prob = plarge)) %>%
    dplyr::mutate(effect_large = rnorm(n = m, mean = 0, sd = sqrt(pge * hsq / sum(large)))) %>%
    dplyr::mutate(effect_small = rnorm(n = m, mean = 0, sd = sqrt((1 - pge) * hsq / (m - sum(large))))) %>%
    dplyr::mutate(effect = large * effect_large + (1 - large) * effect_small) %>%
    dplyr::select(snp_id, effect) %>%
    vroom::vroom_write(file = outfile, col_names = FALSE)
}
# scaled t distribution with 4 df
# The scaled, shifted t distribution has mean mean and variance sd^2 * df/(df-2)
# here, df / (df -2) has value 2
# To get the sum of m independent rvs variances to be hsq, we need an average of hsq / m 
# ie, 2sd^2 = hsq / m; so sd = sqrt(hsq / (2m))
for (hsq in hsq_vec){
  outfile <- paste0("../hapmap3/snp_effects/scenarioIV_scaledt_hsq", hsq, ".txt")
  bim_tib %>%
    dplyr::mutate(large = rbinom(n = m, size = 1, prob = plarge)) %>%
    dplyr::mutate(effect_large = metRology::rt.scaled(n = m, df = 4, mean = 0, sd = sqrt(pge * hsq / (2 * sum(large))))) %>%
    dplyr::mutate(effect_small = metRology::rt.scaled(n = m, df = 4, mean = 0, sd = sqrt((1 - pge) * hsq / (2 * (m - sum(large)))))) %>%
    dplyr::mutate(effect = large * effect_large + (1 - large) * effect_small) %>%
    dplyr::select(snp_id, effect) %>%
    vroom::vroom_write(file = outfile, col_names = FALSE)
}
# laplace distribution
# var = 2 * scale^2
for (hsq in hsq_vec){
  outfile <- paste0("../hapmap3/snp_effects/scenarioIV_laplace_hsq", hsq, ".txt")
  bim_tib %>%
    dplyr::mutate(large = rbinom(n = m, size = 1, prob = plarge)) %>%
    dplyr::mutate(effect_large = VGAM::rlaplace(n = m, location = 0, scale = sqrt(pge * hsq / (2 * sum(large))))) %>%
    dplyr::mutate(effect_small = VGAM::rlaplace(n = m, location = 0, scale = sqrt((1 - pge) * hsq / (2 * (m - sum(large)))))) %>%
    dplyr::mutate(effect = large * effect_large + (1 - large) * effect_small) %>%
    dplyr::select(snp_id, effect) %>%
    vroom::vroom_write(file = outfile, col_names = FALSE)
}

