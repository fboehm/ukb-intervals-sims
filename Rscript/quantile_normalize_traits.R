#!/usr/bin/env Rscript

hsq <- 0.1
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

# read trait tib
trait_tib <- vroom::vroom(trait_file, col_names = FALSE) %>%
  dplyr::select(-8) # drop the last column which contains only NAs
# qn the entire trait_tib
trait_qn <- preprocessCore::normalize.quantiles(x = as.matrix(trait_tib[, 3:7]))
trait_tib_qn <- trait_tib %>%
  dplyr::select(X1, X2) %>%
  dplyr::bind_cols(tibble::as_tibble(trait_qn))

ids <- vroom::vroom(file = "../dat/chr22.fam", col_names = FALSE)

# read files to get ids for training, test, verification, and validation sets
val_ids <- vroom::vroom(file = "../dat/validation-ids.txt", col_names = FALSE)
verif_ids <- vroom::vroom(file = "../dat/verification-ids.txt", col_names = FALSE)
test_ids <- list()
training_ids <- list()
ref_ids <- vroom::vroom(file = "../dat/reference/01_idx.txt", col_names = FALSE)
# subset to get 
trait_tib_validation_qn <- trait_tib_qn %>%
  dplyr::filter(X1 %in% val_ids$X1)
trait_tib_verification_qn <- trait_tib_qn %>%
  dplyr::filter(X1 %in% verif_ids$X1)
trait_tib_reference_qn <- trait_tib_qn %>%
  dplyr::filter(X1 %in% ref_ids$X1)



for (fold in 1:5){
  test_ids[[fold]] <- vroom::vroom(file = paste0("../dat/test-ids-fold", fold, ".txt"), col_names = FALSE)
  training_ids[[fold]] <- vroom::vroom(file = paste0("../dat/training-ids-fold", fold, ".txt"), col_names = FALSE)
}  
# assemble a single tibble with all qn trait values
qn_tib <- trait_tib_qn %>%
  dplyr::left_join(ids, by = c("X1", "X2")) %>%
  dplyr::select(X1, X2, X3, X4, X5, V1:V5)


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

