library(magrittr)
# use the same validation set as for the 136k sims
val_ids <- vroom::vroom("../hapmap3-136k/subjects_for_sims_validation-136k.txt", col_names = FALSE)

scenario_vec <- c("I", "II", "III", "IV")
distribution_vec <- c("laplace", "normal", "scaledt")
hsq_vec <- c(0.1, 0.2, 0.5)
out_dir <- "../dat-quant-small-scale-simulations/validation/"
if (!dir.exists(out_dir)){
    dir.create(out_dir, recursive = TRUE)
}
for (scenario in scenario_vec){
  for (distribution in distribution_vec){
    for (hsq in hsq_vec){
      trait_file <- paste0("../hapmap3-small-scale-simulations/sim_traits_manual/sims_scenario", scenario, "_", distribution, "_hsq", hsq, ".txt")
      trait_tib <- vroom::vroom(trait_file, col_names = FALSE) 
      for (replicate in 1:10){
        # subset to val ids
        trait_tib %>%
          dplyr::right_join(val_ids, by = c("X1", "X2")) %>%
          dplyr::arrange(X1) %>%
          dplyr::select(2 + replicate) %>%
          vroom::vroom_write(file = paste0(out_dir, "scenario", scenario, "_", distribution, "_hsq", hsq, "_replicate", replicate, ".txt"), col_names = FALSE)
      }
    }
  }
}


