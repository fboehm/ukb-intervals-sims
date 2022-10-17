# GOAL: Make multi-trait fam files for use with gemma
# Each replicate gets exactly one column, so ten columns of traits in total
# we have 36 settings, so 36 fam files are needed

n_folds <- 5
library(magrittr)
ids <- vroom::vroom(file = "../hapmap3/plink_files_for_sims/chr1.fam", col_names = FALSE)

scenario_vec <- c("I", "II", "III", "IV")
distribution_vec <- c("laplace", "normal", "scaledt")
hsq_vec <- c(0.1, 0.2, 0.5)
test_ids <- list()
training_ids <- list()
for (fold in 1:n_folds){
  test_ids[[fold]] <- vroom::vroom(file = paste0("../hapmap3/test-ids-fold", fold,".txt" ), col_names = FALSE)
  training_ids[[fold]] <- vroom::vroom(file = paste0("../hapmap3/training-ids-fold", fold,".txt" ), col_names = FALSE)
}

false_to_na <- function(vec){
  vec[!vec] <- NA
  return(vec)
}



for (scenario in scenario_vec){
  for (distribution in distribution_vec){
    for (hsq in hsq_vec){
      trait_file <- paste0("../hapmap3/sim_traits_manual/sims_scenario", scenario, "_", distribution, "_hsq", hsq, ".txt")
      trait_tib <- vroom::vroom(trait_file, col_names = FALSE) 
      # join two tibbles by ids
      all_fam <- ids %>%
        dplyr::select(-6) %>%
        dplyr::left_join(trait_tib, by = c("X1", "X2"))
      # make binary indicators of membership in training set
      tr_indic <- purrr::map(.x = training_ids, .f = function(x){
        false_to_na(all_fam$X1 %in% x$X1)
        }) %>%
        dplyr::bind_cols() %>%
        dplyr::rename_with(.fn = function(x){
          stringr::str_replace_all(x, pattern = "...", replacement = "train")
          })

      training_fam <- all_fam %>%
        dplyr::rename(V1 = X3.y,
                      V2 = X4.y,
                      V3 = X5.y, 
                      V4 = X6, V5 = X7, V6 = X8, V7 = X9, V8 = X10, V9 = X11, V10 = X12,
                      X3 = X3.x, X4 = X4.x, X5 = X5.x) %>%
        dplyr::mutate(dplyr::across(.cols = V1:V10, 
                                    .fns = 
                                      ~ .x * tr_indic
                                    , 
                                    .names = "{.fn}.{.col}")) %>%
        tidyr::unnest(cols = c(`1.V1`, `1.V2`, `1.V3`, `1.V4`, `1.V5`, `1.V6`, `1.V7`, `1.V8`, 
                             `1.V9`, `1.V10`), names_sep = ".") %>%
        dplyr::select(-V1, - V2, - V3, -V4, -V5, -V6, -V7, -V8, -V9, -V10) 
      training_means <- dplyr::summarise(training_fam, dplyr::across(.cols = 6:55, 
                                    .fns = ~ mean(.x, na.rm = TRUE),
                                    .names = "{.fn}.{.col}"))
      training_sds <- dplyr::summarise(training_fam, 
                                       dplyr::across(.cols = 6:55, 
                                                     .fns = ~ sd(.x, na.rm = TRUE),
                                                     .names = "{.fn}.{.col}"))
# 
      te_indic <- purrr::map(.x = test_ids, .f = function(x){false_to_na(all_fam$X1 %in% x$X1)}) %>%
        dplyr::bind_cols() %>%
        dplyr::rename_with(.fn = function(x){
          stringr::str_replace_all(x, pattern = "...", replacement = "test")
        })
      test_fam <- all_fam %>%
        dplyr::rename(V1 = X3.y,
                      V2 = X4.y,
                      V3 = X5.y, 
                      V4 = X6, V5 = X7, V6 = X8, V7 = X9, V8 = X10, V9 = X11, V10 = X12,
                      X3 = X3.x, X4 = X4.x, X5 = X5.x) %>%
        dplyr::mutate(dplyr::across(.cols = V1:V10, 
                                    .fns = 
                                      ~ .x * te_indic
                                    , 
                                    .names = "{.fn}.{.col}")) %>%
        tidyr::unnest(cols = c(`1.V1`, `1.V2`, `1.V3`, `1.V4`, `1.V5`, `1.V6`, `1.V7`, `1.V8`, 
                               `1.V9`, `1.V10`), names_sep = ".") %>%
        dplyr::select(-V1, - V2, - V3, -V4, -V5, -V6, -V7, -V8, -V9, -V10)
      # center & scale the test fam traits
      te_centered <- test_fam %>%
        dplyr::select(- c(1:5)) %>%
        as.matrix() %>%
        sweep(MARGIN = 2, STATS = as.array(unlist(training_means)), FUN = "-")
      te_scaled <- te_centered %>%
        sweep(MARGIN = 2, STATS = as.array(unlist(training_sds)), FUN = "/")
      test_fam2 <- test_fam %>%
        dplyr::select(1:5) %>%
        dplyr::bind_cols(te_scaled)
      # center & scale the training fam traits
      tr_centered <- training_fam %>%
        dplyr::select(- c(1:5)) %>%
        as.matrix() %>%
        sweep(MARGIN = 2, STATS = as.array(unlist(training_means)), FUN = "-")
      tr_scaled <- tr_centered %>%
        sweep(MARGIN = 2, STATS = as.array(unlist(training_sds)), FUN = "/")
      training_fam2 <- training_fam %>%
        dplyr::select(1:5) %>%
        dplyr::bind_cols(tr_scaled)
      
      
      
      # write resulting fam file to appropriate subdir
      outdir <- paste0("../dat-quant-5fold-manual-sim/gemma/scenario", scenario, "/", distribution, "/hsq", hsq)
      if (!dir.exists(outdir)){
        dir.create(outdir, recursive = TRUE)
      }
      out_fn <- paste0(outdir, "/chr1.fam")
      training_fam2 %>%
        vroom::vroom_write(file = out_fn, col_names = FALSE)
      test_out_fn <- paste0(outdir, "/test.fam")
      test_fam2 %>%
        vroom::vroom_write(file = test_out_fn, col_names = FALSE)
    }
  }
}
