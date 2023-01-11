# Goal here is to reuse the verification set subjects from the sims with 136,000 subjects
# to examine the interval widths.

library(magrittr)

ids_file <- "../hapmap3-136k/subjects_for_sims-136k.txt"
sim_ids <- vroom::vroom(file = ids_file, col_names = FALSE)
# sample verification 1000 subjects
# write verif_ids to file
verif_ids_file <- "../hapmap3-136k/subjects_for_sims_verification-136k.txt"
verif_ids <- vroom::vroom(file = verif_ids_file, col_names = FALSE)
# take setdiff to get 135,000 remaining ids
remaining_ids <- setdiff(sim_ids$X1, verif_ids$X1)
# read 10000 validation set
valid_ids_file <- "../hapmap3-136k/subjects_for_sims_validation-136k.txt"
valid_ids <- vroom::vroom(file = valid_ids_file, col_names = FALSE)
# get test&training 125000 subjects
test_and_training_ids_big <- setdiff(remaining_ids, valid_ids$X1)
# sample from test_and_training_big to get 12,500 subjects
set.seed(2023-01-11)
test_and_training_ids <- sample(x = test_and_training_ids_big, size = 12500, replace = FALSE)
# set up for 5-fold cv with the remaining 11,000 subjects
ids_shuffled <- sample(x = test_and_training_ids, size = length(test_and_training_ids))
folds <- cut(seq(1,length(ids_shuffled)),breaks=5,labels=FALSE)
# https://stats.stackexchange.com/questions/61090/how-to-split-a-data-set-to-do-10-fold-cross-validation
for (i in 1:5){
  ids_tib <- tibble::tibble(ids_shuffled, folds) %>%
    dplyr::arrange(ids_shuffled)
  # write test ids for the fold
  ids_tib %>%
    dplyr::filter(folds == i) %>%
    dplyr::mutate(X2 = ids_shuffled) %>%
    dplyr::select(- folds) %>%
    #vroom::vroom_write(file = snakemake@output[[2 + i]], col_names = FALSE)
    vroom::vroom_write(file = paste0("../hapmap3-small-scale-simulations/test-ids-fold", i, ".txt"), col_names = FALSE)
  # write training ids for the fold
  ids_tib %>%
    dplyr::filter(folds != i) %>%
    dplyr::mutate(X2 = ids_shuffled) %>%
    dplyr::select(- folds) %>%
    vroom::vroom_write(file = paste0("../hapmap3-small-scale-simulations/training-ids-fold", i, ".txt"), col_names = FALSE)
}
# lastly, combine all subjects into a single file for use with plink.
all_ids <- c(test_and_training_ids, verif_ids$X1, valid_ids$X1)
out_file <- "../hapmap3-small-scale-simulations/subjects_for_sims.txt"
tibble::tibble(X1 = all_ids, X2 = all_ids) %>%
    dplyr::arrange(X1) %>%
    vroom::vroom_write(file = out_file, col_names = FALSE)



