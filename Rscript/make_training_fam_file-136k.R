# GOAL: Make multi-trait fam files for use with gemma
# Each replicate gets exactly one column, so ten columns of traits in total
# we have 36 settings, so 36 fam files are needed

library(magrittr)
ids <- vroom::vroom(file = "../hapmap3-136k/plink_files_for_sims/chr1.fam", col_names = FALSE)

scenario_vec <- c("I", "II", "III", "IV")
distribution_vec <- c("laplace", "normal", "scaledt")
hsq_vec <- c(0.1, 0.2, 0.5)
test_ids <- list()
training_ids <- list()
for (fold in 1:5){
  test_ids[[fold]] <- vroom::vroom(file = paste0("../hapmap3-136k/test-ids-fold", fold,"-136k.txt" ), col_names = FALSE)
  training_ids[[fold]] <- vroom::vroom(file = paste0("../hapmap3-136k/training-ids-fold", fold,"-136k.txt" ), col_names = FALSE)
}

false_to_na <- function(vec){
  vec[!vec] <- NA
  return(vec)
}



for (scenario in scenario_vec){
  for (distribution in distribution_vec){
    for (hsq in hsq_vec){
      trait_file <- paste0("../hapmap3-136k/sim_traits/sims_scenario", scenario, "_", distribution, "_hsq", hsq, ".txt.phen")
      trait_tib <- vroom::vroom(trait_file, col_names = FALSE) %>%
        dplyr::select(-13) # drop the last column which contains only NAs
      # join two tibbles by ids
      all_fam <- ids %>%
        dplyr::select(-6) %>%
        dplyr::left_join(trait_tib, by = c("X1", "X2"))
      
      
      
      
      
      # make binary indicators of membership in training set
      # make binary indicators of membership in training set
      tr_indic <- tibble::tibble(train1 = all_fam$X1 %in% training_ids[[1]]$X1,
                                 train2 = all_fam$X1 %in% training_ids[[2]]$X1,
                                 train3 = all_fam$X1 %in% training_ids[[3]]$X1,
                                 train4 = all_fam$X1 %in% training_ids[[4]]$X1,
                                 train5 = all_fam$X1 %in% training_ids[[5]]$X1
      ) %>%
        dplyr::mutate(fold1_na = false_to_na(train1),
                      fold2_na = false_to_na(train2),
                      fold3_na = false_to_na(train3),
                      fold4_na = false_to_na(train4),
                      fold5_na = false_to_na(train5)
        )
      training_fam <- all_fam %>%
        dplyr::rename(V1 = X3.y,
                      V2 = X4.y,
                      V3 = X5.y, 
                      V4 = X6, V5 = X7, V6 = X8, V7 = X9, V8 = X10, V9 = X11, V10 = X12,
                      X3 = X3.x, X4 = X4.x, X5 = X5.x) %>%
        dplyr::mutate(tr1_fold1 = V1 * tr_indic$fold1_na,
                      tr1_fold2 = V1 * tr_indic$fold2_na,
                      tr1_fold3 = V1 * tr_indic$fold3_na,
                      tr1_fold4 = V1 * tr_indic$fold4_na,
                      tr1_fold5 = V1 * tr_indic$fold5_na,
                      tr2_fold1 = V2 * tr_indic$fold1_na,
                      tr2_fold2 = V2 * tr_indic$fold2_na,
                      tr2_fold3 = V2 * tr_indic$fold3_na,
                      tr2_fold4 = V2 * tr_indic$fold4_na,
                      tr2_fold5 = V2 * tr_indic$fold5_na,
                      tr3_fold1 = V3 * tr_indic$fold1_na,
                      tr3_fold2 = V3 * tr_indic$fold2_na,
                      tr3_fold3 = V3 * tr_indic$fold3_na,
                      tr3_fold4 = V3 * tr_indic$fold4_na,
                      tr3_fold5 = V3 * tr_indic$fold5_na,
                      tr4_fold1 = V4 * tr_indic$fold1_na,
                      tr4_fold2 = V4 * tr_indic$fold2_na,
                      tr4_fold3 = V4 * tr_indic$fold3_na,
                      tr4_fold4 = V4 * tr_indic$fold4_na,
                      tr4_fold5 = V4 * tr_indic$fold5_na,
                      tr5_fold1 = V5 * tr_indic$fold1_na,
                      tr5_fold2 = V5 * tr_indic$fold2_na,
                      tr5_fold3 = V5 * tr_indic$fold3_na,
                      tr5_fold4 = V5 * tr_indic$fold4_na,
                      tr5_fold5 = V5 * tr_indic$fold5_na,
                      tr6_fold1 = V6 * tr_indic$fold1_na,
                      tr6_fold2 = V6 * tr_indic$fold2_na,
                      tr6_fold3 = V6 * tr_indic$fold3_na,
                      tr6_fold4 = V6 * tr_indic$fold4_na,
                      tr6_fold5 = V6 * tr_indic$fold5_na,
                      tr7_fold1 = V7 * tr_indic$fold1_na,
                      tr7_fold2 = V7 * tr_indic$fold2_na,
                      tr7_fold3 = V7 * tr_indic$fold3_na,
                      tr7_fold4 = V7 * tr_indic$fold4_na,
                      tr7_fold5 = V7 * tr_indic$fold5_na,
                      tr8_fold1 = V8 * tr_indic$fold1_na,
                      tr8_fold2 = V8 * tr_indic$fold2_na,
                      tr8_fold3 = V8 * tr_indic$fold3_na,
                      tr8_fold4 = V8 * tr_indic$fold4_na,
                      tr8_fold5 = V8 * tr_indic$fold5_na,
                      tr9_fold1 = V9 * tr_indic$fold1_na,
                      tr9_fold2 = V9 * tr_indic$fold2_na,
                      tr9_fold3 = V9 * tr_indic$fold3_na,
                      tr9_fold4 = V9 * tr_indic$fold4_na,
                      tr9_fold5 = V9 * tr_indic$fold5_na,
                      tr10_fold1 = V10 * tr_indic$fold1_na,
                      tr10_fold2 = V10 * tr_indic$fold2_na,
                      tr10_fold3 = V10 * tr_indic$fold3_na,
                      tr10_fold4 = V10 * tr_indic$fold4_na,
                      tr10_fold5 = V10 * tr_indic$fold5_na,
        ) %>%
        dplyr::select(-V1, - V2, - V3, -V4, -V5, -V6, -V7, -V8, -V9, -V10) 
      
      # write resulting fam file to appropriate subdir
      outdir <- paste0("../dat-quant-136k/gemma/scenario", scenario, "/", distribution, "/hsq", hsq)
      if (!dir.exists(outdir)){
        dir.create(outdir, recursive = TRUE)
      }
      out_fn <- paste0(outdir, "/chr1.fam")
      training_fam %>%
        vroom::vroom_write(file = out_fn, col_names = FALSE)
      
    }
  }
}
