#!/usr/bin/env Rscript

hsq <- 0.5
pcausal <- 0.1

library(magrittr)

false_to_na <- function(vec){
  vec[!vec] <- NA
  return(vec)
}

trait_file <- paste0("../dat/hsq", hsq, "_pcausal", pcausal, "/sim_traits/sims_Chr22_hsq", hsq, "_pcausal", pcausal, ".phen")
output_with_nas <- paste0("../dat/hsq", hsq, "_pcausal", pcausal, "/sim_traits/sims_Chr22_hsq", hsq, "_pcausal", pcausal, "-NAs.fam")
output_without_nas <- paste0("../dat/hsq", hsq, "_pcausal", pcausal, "/sim_traits/fam_all_subjects/qn_traits_all_subjects_hsq", hsq, "_pcausal", pcausal, ".fam")

# make directories needed
output_dir <- paste0("../dat/hsq", hsq, "_pcausal", pcausal, "/sim_traits/fam_all_subjects")
if (!dir.exists(output_dir)){
  dir.create(output_dir)
} else {
  print("Dir already exists!")
}
output_dir <- paste0("../dat/hsq", hsq, "_pcausal", pcausal, "/validation")
if (!dir.exists(output_dir)){
  dir.create(output_dir)
} else {
  print("Dir already exists!")
}

output_dir <- paste0("../dat/hsq", hsq, "_pcausal", pcausal, "/verification")
if (!dir.exists(output_dir)){
  dir.create(output_dir)
} else {
  print("Dir already exists!")
}

                             
trait_tib <- vroom::vroom(trait_file, col_names = FALSE) %>%
  dplyr::select(-8) # drop the last column which contains only NAs

ids <- vroom::vroom(file = "../dat/chr22.fam", col_names = FALSE)


# read files to get ids for training, test, and validation sets
val_ids <- vroom::vroom(file = paste0("../dat/hsq", hsq, "_pcausal", pcausal, "/validation-ids.txt"), col_names = FALSE)
verif_ids <- vroom::vroom(file = paste0("../dat/hsq", hsq, "_pcausal", pcausal, "/verification-ids.txt"), col_names = FALSE)
test_ids <- list()
training_ids <- list()
#
for (fold in 1:5){
  test_ids[[fold]] <- vroom::vroom(file = paste0("../dat/hsq", hsq, "_pcausal", pcausal, "/test-ids-fold", fold, ".txt"), col_names = FALSE)
  training_ids[[fold]] <- vroom::vroom(file = paste0("../dat/hsq", hsq, "_pcausal", pcausal, "/training-ids-fold", fold, ".txt"), col_names = FALSE)
  
  # use training set only to determine the distribution for quantile normalization
  # then, use the inferred distribution to quantile normalize, separately, the training and the validation and the verification and the test sets.
  trait_tib_training <- trait_tib %>%
    dplyr::filter(X1 %in% training_ids[[fold]]$X1)
  trait_tib_test <- trait_tib %>%
    dplyr::filter(X1 %in% test_ids[[fold]]$X1)
  trait_tib_validation <- trait_tib %>%
    dplyr::filter(X1 %in% val_ids$X1)
  trait_tib_verification <- trait_tib %>%
    dplyr::filter(X1 %in% verif_ids$X1)
  
  
  training_qn <- preprocessCore::normalize.quantiles(as.matrix(trait_tib_training[, 3:7]))
  target_test <- preprocessCore::normalize.quantiles.determine.target(x = training_qn, target.length = nrow(trait_tib_test))
  test_qn <- preprocessCore::normalize.quantiles.use.target(x = as.matrix(trait_tib_test[, 3:7]), target = target_test)
  target_validation <- preprocessCore::normalize.quantiles.determine.target(x = training_qn, target.length = nrow(trait_tib_validation))
  validation_qn <- preprocessCore::normalize.quantiles.use.target(x = as.matrix(trait_tib_validation[, 3:7]), target = target_validation)
  target_verification <- preprocessCore::normalize.quantiles.determine.target(x = training_qn, target.length = nrow(trait_tib_verification))
  verification_qn <- preprocessCore::normalize.quantiles.use.target(x = as.matrix(trait_tib_verification[, 3:7]), target = target_verification)
  
  ######
  training_qn_tib <- training_ids[[fold]] %>%
    dplyr::bind_cols(tibble::as_tibble(training_qn))
  test_qn_tib <- test_ids[[fold]] %>%
    dplyr::bind_cols(tibble::as_tibble(test_qn))
  validation_qn_tib <- val_ids %>%
    dplyr::bind_cols(tibble::as_tibble(validation_qn))
  verification_qn_tib <- verif_ids %>%
    dplyr::bind_cols(tibble::as_tibble(verification_qn))
  
  qn_tib_pre <- training_qn_tib %>%
    dplyr::bind_rows(test_qn_tib) %>%
    dplyr::bind_rows(validation_qn_tib) %>%
    dplyr::bind_rows(verification_qn_tib) %>%
    dplyr::arrange(X1) #%>%
  qn_tib <- qn_tib_pre %>%
    dplyr::left_join(ids, by = c("X1", "X2")) %>%
    dplyr::select(X1, X2, X3, X4, X5, V1:V5)
} # end loop over folds


# make binary indicators of membership in training set
tr_indic <- tibble::tibble(train1 = qn_tib$X1 %in% training_ids[[1]]$X1,
                           train2 = qn_tib$X1 %in% training_ids[[2]]$X1,
                           train3 = qn_tib$X1 %in% training_ids[[3]]$X1,
                           train4 = qn_tib$X1 %in% training_ids[[4]]$X1,
                           train5 = qn_tib$X1 %in% training_ids[[5]]$X1
) %>%
  dplyr::mutate(fold1_na = false_to_na(train1),
                fold2_na = false_to_na(train2),
                fold3_na = false_to_na(train3),
                fold4_na = false_to_na(train4),
                fold5_na = false_to_na(train5)
  )
qn_tib2 <- qn_tib %>%
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
                tr5_fold5 = V5 * tr_indic$fold5_na
  ) %>%
  dplyr::select(-V1, - V2, - V3, -V4, -V5) 

#####
qn_tib2 %>%
  vroom::vroom_write(file = output_with_nas, col_names = FALSE)
# write qn traits for all subjects in a single fam file
qn_tib %>%
  vroom::vroom_write(file = output_without_nas, col_names = FALSE)

# make the pheno files for verification set
qt_verif <- qn_tib %>%
  dplyr::arrange(X1) %>%
  dplyr::filter(X1 %in% verif_ids$X1)
for (col_num in 6:10){
  outfile <- paste0("../dat/hsq", hsq, "_pcausal", pcausal, "/verification/pheno", col_num - 5, "_hsq", hsq, "_pcausal", pcausal, ".txt")
  qt_verif %>%
    dplyr::select(dplyr::all_of(col_num)) %>%
    vroom::vroom_write(file = outfile, col_names = FALSE)
}

# make the pheno files for validation set
qt_valid <- qn_tib %>%
  dplyr::arrange(X1) %>%
  dplyr::filter(X1 %in% val_ids$X1)
for (col_num in 6:10){
  outfile <- paste0("../dat/hsq", hsq, "_pcausal", pcausal, "/validation/pheno", col_num - 5, "_hsq", hsq, "_pcausal", pcausal, ".txt")
  qt_valid %>%
    dplyr::select(dplyr::all_of(col_num)) %>%
    vroom::vroom_write(file = outfile, col_names = FALSE)
}

