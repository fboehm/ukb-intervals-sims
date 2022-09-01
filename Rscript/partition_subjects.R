#hsq <- 0.5
#pcausal <- 0.1

library(magrittr)

#ids <- vroom::vroom(file = snakemake@input[[1]], col_names = FALSE)
## MAKE REFERENCE DATA SUBSET
load("/net/mulan/disk2/yasheng/comparisonProject/02_pheno/01_sqc.RData")
ref_num <- 250 # reference number: 250 males, 250 females 
# male
idx_male <- which(sqc_i$Inferred.Gender == "M")
sqc_male <- sqc_i[idx_male, c("Inferred.Gender", "idx")]
set.seed(20170529)
idx_male_ref <- sample(sqc_male$idx, ref_num)
# female
idx_female <- which(sqc_i$Inferred.Gender == "F")
sqc_female <- sqc_i[idx_female, c("Inferred.Gender", "idx")]
set.seed(20170529)
idx_female_ref <- sample(sqc_female$idx, ref_num)
# index for reference 
idx_ref <- sort(c(idx_male_ref, idx_female_ref))
write.table(cbind(idx_ref, idx_ref), 
             file = "../dat/reference/01_idx.txt", 
             col.names = FALSE, row.names = FALSE, quote = FALSE)
sqc_sub <- sqc_i[!sqc_i$idx %in% idx_ref,]
##
set.seed(2022-06-01)
pre_val_ids <- sample(x = sqc_sub$idx, size = 36629, replace = FALSE) 
# split pre_val_ids into two - one validation set and a "verification" set of, say, 1000.
# verification set is used in the cv, as the X_{n+1}, while validation is used in the tuning of parameters for DBSLMM
# we also need the "reference" set
verif_ids <- sample(pre_val_ids, size = 1000, replace = FALSE)
val_ids <- setdiff(pre_val_ids, verif_ids)
# Write to files
val_ids %>%
  tibble::tibble(.name_repair = "universal") %>%
  dplyr::rename(X1 = 1) %>%
  dplyr::arrange(X1) %>%
  dplyr::mutate(X2 = X1) %>%
  #vroom::vroom_write(file = snakemake@output[[1]], col_names = FALSE)
  vroom::vroom_write(file = "../dat/validation-ids.txt", col_names = FALSE)
verif_ids %>%
  tibble::tibble(.name_repair = "universal") %>%
  dplyr::rename(X1 = 1) %>%
  dplyr::arrange(X1) %>%
  dplyr::mutate(X2 = X1) %>%
  #vroom::vroom_write(file = snakemake@output[[2]], col_names = FALSE)
  vroom::vroom_write(file = "../dat/verification-ids.txt", col_names = FALSE)

test_and_training_ids <- setdiff(sqc_sub$idx, pre_val_ids)
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
    vroom::vroom_write(file = paste0("../dat/test-ids-fold", i, ".txt"), col_names = FALSE)
  # write training ids for the fold
  ids_tib %>%
    dplyr::filter(folds != i) %>%
    dplyr::mutate(X2 = ids_shuffled) %>%
    dplyr::select(- folds) %>%
    vroom::vroom_write(file = paste0("../dat/training-ids-fold", i, ".txt"), col_names = FALSE)
}

## Verify subsetting
intersect(test_and_training_ids, idx_ref)
intersect(test_and_training_ids, val_ids)
intersect(test_and_training_ids, verif_ids)
intersect(verif_ids, val_ids)
intersect(verif_ids, idx_ref)
intersect(val_ids, idx_ref)
