# GOAL: Make multi-trait fam files for use with gemma
# Each replicate gets exactly one column, so ten columns of traits in total
# we have 36 settings, so 36 fam files are needed

library(magrittr)
ids <- vroom::vroom(file = "../hapmap3/plink_files_for_sims/chr1.fam", col_names = FALSE)

scenario_vec <- c("I", "II", "III", "IV")
distribution_vec <- c("laplace", "normal", "scaledt")
hsq_vec <- c(0.1, 0.2, 0.5)
test_ids <- list()
training_ids <- list()
for (rep in 1:10){
  test_ids[[rep]] <- vroom::vroom(file = paste0("../hapmap3/subjects_for_sims_test_replicate", rep, ".txt"), col_names = FALSE)
  training_ids[[rep]] <- vroom::vroom(file = paste0("../hapmap3/subjects_for_sims_training_replicate", rep, ".txt"), col_names = FALSE)
}  

false_to_na <- function(vec){
  vec[!vec] <- NA
  return(vec)
}



for (scenario in scenario_vec){
  for (distribution in distribution_vec){
    for (hsq in hsq_vec){
      trait_file <- paste0("../hapmap3/sim_traits/sims_scenario", scenario, "_", distribution, "_hsq", hsq, ".txt.phen")
      trait_tib <- vroom::vroom(trait_file, col_names = FALSE) %>%
        dplyr::select(-13) # drop the last column which contains only NAs
      # join two tibbles by ids
      all_fam <- ids %>%
        dplyr::select(-6) %>%
        dplyr::left_join(trait_tib, by = c("X1", "X2"))
      
      
      
      
      
      # make binary indicators of membership in training set
      tr_indic <- tibble::tibble(train1 = all_fam$X1 %in% training_ids[[1]]$X1,
                                 train2 = all_fam$X1 %in% training_ids[[2]]$X1,
                                 train3 = all_fam$X1 %in% training_ids[[3]]$X1,
                                 train4 = all_fam$X1 %in% training_ids[[4]]$X1,
                                 train5 = all_fam$X1 %in% training_ids[[5]]$X1,
                                 train6 = all_fam$X1 %in% training_ids[[6]]$X1,
                                 train7 = all_fam$X1 %in% training_ids[[7]]$X1,
                                 train8 = all_fam$X1 %in% training_ids[[8]]$X1,
                                 train9 = all_fam$X1 %in% training_ids[[9]]$X1,
                                 train10 = all_fam$X1 %in% training_ids[[10]]$X1
      ) %>%
        dplyr::mutate(rep1_na = false_to_na(train1),
                      rep2_na = false_to_na(train2),
                      rep3_na = false_to_na(train3),
                      rep4_na = false_to_na(train4),
                      rep5_na = false_to_na(train5),
                      rep6_na = false_to_na(train6),
                      rep7_na = false_to_na(train7),
                      rep8_na = false_to_na(train8),
                      rep9_na = false_to_na(train9),
                      rep10_na = false_to_na(train10)
        )
      training_fam <- all_fam %>%
        dplyr::mutate(rep1 = X3.y * tr_indic$rep1_na,
                      rep2 = X4.y * tr_indic$rep2_na,
                      rep3 = X5.y * tr_indic$rep3_na,
                      rep4 = X6 * tr_indic$rep4_na,
                      rep5 = X7 * tr_indic$rep5_na,
                      rep6 = X8 * tr_indic$rep6_na,
                      rep7 = X9 * tr_indic$rep7_na,
                      rep8 = X10 * tr_indic$rep8_na,
                      rep9 = X11 * tr_indic$rep9_na,
                      rep10 = X12 * tr_indic$rep10_na
        ) %>%
        dplyr::select(-(6:15))
      # write resulting fam file to appropriate subdir
      out_fn <- paste0("../dat-quant/gemma/scenario", scenario, "/", distribution, "/hsq", hsq, "/chr1.fam")
      training_fam %>%
        vroom::vroom_write(file = out_fn)
      
    }
  }
}
