# GOAL: choose 13,000 subjects for simulating traits with.
# These 13,000 include 1000 for validation, 1000 for verification, and 
# 11,000 for training and test together
# we'll do five-fold crossvalidation with the 11,000 - 1000 test and 10,000 training each time


library(magrittr)
# sample 13,000 subjects from the 337,129
## read fam file for entire 337,129
fam_file <- "../hapmap3/hm3_ukb/chr1.fam"
fam_tib <- vroom::vroom(fam_file, col_names = FALSE)
set.seed(2022-09-14)
# sample 13000 ids
sim_ids <- sample(x = fam_tib$X1, size = 13000, replace = FALSE)
# write 13000 ids for use with plink --keep
ids_file <- "../hapmap3/subjects_for_sims.txt"
tibble::tibble(X1 = sim_ids, X2 = sim_ids) %>%
  dplyr::arrange(X1) %>%
  vroom::vroom_write(file = ids_file, col_names = FALSE)
# sample verification 1000 subjects
verif_ids <- sample(x = sim_ids, size = 1000, replace = FALSE)
# write verif_ids to file
verif_ids_file <- "../hapmap3/subjects_for_sims_verification.txt"
tibble::tibble(X1 = verif_ids, X2 = verif_ids) %>%
  dplyr::arrange(X1) %>%
  vroom::vroom_write(file = verif_ids_file, col_names = FALSE)
# take setdiff to get 12,000 remaining ids
remaining_ids <- setdiff(sim_ids, verif_ids)
# sample 1000 validation set
valid_ids <- sample(x = remaining_ids, size = 1000, replace = FALSE)
# write validation ids to file
valid_ids_file <- "../hapmap3/subjects_for_sims_validation.txt"
tibble::tibble(X1 = valid_ids, X2 = valid_ids) %>%
  dplyr::arrange(X1) %>%
  vroom::vroom_write(file = valid_ids_file, col_names = FALSE)
# get test&training 11000 subjects
test_and_training_ids <- setdiff(remaining_ids, valid_ids)

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
    vroom::vroom_write(file = paste0("../hapmap3/test-ids-fold", i, ".txt"), col_names = FALSE)
  # write training ids for the fold
  ids_tib %>%
    dplyr::filter(folds != i) %>%
    dplyr::mutate(X2 = ids_shuffled) %>%
    dplyr::select(- folds) %>%
    vroom::vroom_write(file = paste0("../hapmap3/training-ids-fold", i, ".txt"), col_names = FALSE)
}
# lastly, choose 500 of the validation subjects to be the reference panel
ref_ids <- sample(x = valid_ids, size = 500, replace = FALSE)
ref_ids_file <- "../hapmap3/subjects_for_sims_reference.txt"
tibble::tibble(X1 = ref_ids, X2 = ref_ids) %>%
  dplyr::arrange(X1) %>%
  vroom::vroom_write(file = ref_ids_file, col_names = FALSE)



