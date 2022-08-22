hsq <- 0.5
pcausal <- 0.1

library(magrittr)

set.seed(2022-06-01)
#ids <- vroom::vroom(file = snakemake@input[[1]], col_names = FALSE)
ids <- vroom::vroom(file = "../dat/chr22.fam", col_names = FALSE)
pre_val_ids <- sample(x = ids$X1, size = 37129, replace = FALSE) 
# split pre_val_ids into two - one validation set and a "verification" set of, say, 1000.
# verification set is used in the cv, as the X_{n+1}, while validation is used in the tuning of parameters for DBSLMM
verif_ids <- sample(pre_val_ids, size = 1000, replace = FALSE)
val_ids <- setdiff(pre_val_ids, verif_ids)
# Write to files
val_ids %>%
  tibble::tibble(.name_repair = "universal") %>%
  dplyr::rename(X1 = 1) %>%
  dplyr::arrange(X1) %>%
  dplyr::mutate(X2 = X1) %>%
  #vroom::vroom_write(file = snakemake@output[[1]], col_names = FALSE)
  vroom::vroom_write(file = paste0("../dat/hsq", hsq, "_pcausal", pcausal, "/validation-ids.txt"), col_names = FALSE)
verif_ids %>%
  tibble::tibble(.name_repair = "universal") %>%
  dplyr::rename(X1 = 1) %>%
  dplyr::arrange(X1) %>%
  dplyr::mutate(X2 = X1) %>%
  #vroom::vroom_write(file = snakemake@output[[2]], col_names = FALSE)
  vroom::vroom_write(file = paste0("../dat/hsq", hsq, "_pcausal", pcausal, "/verification-ids.txt"), col_names = FALSE)

test_and_training_ids <- setdiff(ids$X1, pre_val_ids)
# set up for 5-fold cv with the remaining 300,000 subjects
set.seed(2022-06-10)
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
    vroom::vroom_write(file = paste0("../dat/hsq", hsq, "_pcausal", pcausal, "/test-ids-fold", i, ".txt"), col_names = FALSE)
  # write training ids for the fold
  ids_tib %>%
    dplyr::filter(folds != i) %>%
    dplyr::mutate(X2 = ids_shuffled) %>%
    dplyr::select(- folds) %>%
    vroom::vroom_write(file = paste0("../dat/hsq", hsq, "_pcausal", pcausal, "/training-ids-fold", i, ".txt"), col_names = FALSE)
}
