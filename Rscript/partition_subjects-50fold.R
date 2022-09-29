# GOAL: choose 13,000 subjects for simulating traits with.
# These 14,500 include 1000 for validation, 1000 for verification, and 
# we'll do fifty-fold crossvalidation with the 12,500 - 250 test and 12,250 training each time


library(magrittr)
# sample 14,500 subjects from the 337,129
## read fam file for entire 337,129
fam_file <- "../hapmap3/hm3_ukb/chr1.fam"
fam_tib <- vroom::vroom(fam_file, col_names = FALSE)
set.seed(2022-09-14)
# sample 13000 ids
#sim_ids <- sample(x = fam_tib$X1, size = 14500, replace = FALSE)
# write 14500 ids for use with plink --keep
ids_file <- "../hapmap3/subjects_for_sims.txt"
sim_ids <- vroom::vroom(ids_file, col_names = FALSE)

# sample verification 1000 subjects
#verif_ids <- sample(x = sim_ids, size = 1000, replace = FALSE)
# write verif_ids to file
verif_ids_file <- "../hapmap3/subjects_for_sims_verification.txt"
verif_ids <- vroom::vroom(verif_ids_file, col_names = FALSE)
valid_ids_file <- "../hapmap3/subjects_for_sims_validation.txt"
valid_ids <- vroom::vroom(valid_ids_file, col_names = FALSE)
# get test&training 12500 subjects

test_and_training_ids <- setdiff(sim_ids, union(valid_ids, verif_ids))
n_folds <- 50
# set up for 50-fold cv with the remaining 11,000 subjects
set.seed(2022-09-29)
ids_shuffled <- sample(x = test_and_training_ids$X1, size = length(test_and_training_ids$X1))
folds <- cut(seq(1,length(ids_shuffled)),breaks=n_folds,labels=FALSE)
# https://stats.stackexchange.com/questions/61090/how-to-split-a-data-set-to-do-10-fold-cross-validation
for (i in 1:n_folds){
  ids_tib <- tibble::tibble(ids_shuffled, folds) %>%
    dplyr::arrange(ids_shuffled)
  # write test ids for the fold
  ids_tib %>%
    dplyr::filter(folds == i) %>%
    dplyr::mutate(X2 = ids_shuffled) %>%
    dplyr::select(- folds) %>%
    #vroom::vroom_write(file = snakemake@output[[2 + i]], col_names = FALSE)
    vroom::vroom_write(file = paste0("../hapmap3-50fold/test-ids-fold", i, ".txt"), col_names = FALSE)
  # write training ids for the fold
  ids_tib %>%
    dplyr::filter(folds != i) %>%
    dplyr::mutate(X2 = ids_shuffled) %>%
    dplyr::select(- folds) %>%
    vroom::vroom_write(file = paste0("../hapmap3-50fold/training-ids-fold", i, ".txt"), col_names = FALSE)
}

