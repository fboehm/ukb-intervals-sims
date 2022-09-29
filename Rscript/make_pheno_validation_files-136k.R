library(magrittr)
val_ids <- vroom::vroom("../hapmap3-136k/subjects_for_sims_validation-136k.txt", col_names = FALSE)

scenario_vec <- c("I", "II", "III", "IV")
distribution_vec <- c("laplace", "normal", "scaledt")
hsq_vec <- c(0.1, 0.2, 0.5)
for (scenario in scenario_vec){
  for (distribution in distribution_vec){
    for (hsq in hsq_vec){
      trait_file <- paste0("../hapmap3-136k/sim_traits/sims_scenario", scenario, "_", distribution, "_hsq", hsq, ".txt.phen")
      trait_tib <- vroom::vroom(trait_file, col_names = FALSE) %>%
        dplyr::select(-13) # drop the last column which contains only NAs
      for (replicate in 1:10){
        # subset to val ids
        trait_tib %>%
          dplyr::right_join(val_ids, by = c("X1", "X2")) %>%
          dplyr::arrange(X1) %>%
          dplyr::select(2 + replicate) %>%
          vroom::vroom_write(file = paste0("../dat-quant-136k/validation/scenario", scenario, "_", distribution, "_hsq", hsq, "_replicate", replicate, ".txt"), col_names = FALSE)
      }
    }
  }
}


