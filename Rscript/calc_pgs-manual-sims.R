# Goal: use DBSLMM-estimated SNP effects and genotypes to calculate PGS
args = commandArgs(trailingOnly=TRUE)
(scenario <- args[1])
(distribution <- args[2])
(hsq <- as.numeric(args[3]))
(pheno <- as.numeric(args[4]))
library(magrittr)

n_folds <- 5
alpha <- 0.2
b_file <- "../hapmap3/plink_files_for_sims/chr1.rds"
gg <- bigsnpr::snp_attach(b_file)
gg_imputed <- bigsnpr::snp_fastImputeSimple(Gna = gg$genotypes) 
# read ids for verification set
v_file <- "../hapmap3/subjects_for_sims_verification.txt"
verif_ids <- readr::read_table(v_file, col_names = FALSE)

# read phen file with "true" phenotypes
test_fam_fn <- paste0("../dat-quant-5fold-manual-sim/gemma/scenario", scenario, "/", distribution, "/hsq", hsq, "/test.fam")
test_fam <- vroom::vroom(test_fam_fn, col_names = FALSE)
coverage <- numeric() 
pgs_list <- list()
muhat <- list() # to store the muhat (Xn+1) values from Eqn 3.1
for (fold in 1:n_folds){
  # DBSLMM outputted file
  d_file <- paste0("../dat-quant-5fold-manual-sim/DBSLMM/scenario", scenario, "/", distribution, "/hsq", hsq, "/summary_ukb_pheno", pheno, "_scenario", scenario, "_", distribution, "_hsq", hsq, "_fold", fold, "_chr1_best.dbslmm.txt")
  dbslmm_output <- readr::read_table(d_file, col_names = FALSE)
  # read plink genotypes bed file
  # read training set membership file
  training_membership_file <- paste0("../hapmap3/training-ids-fold", fold, ".txt")
  training_membership <- readr::read_table(training_membership_file, col_names = FALSE)
  tr_indic <- gg$fam$family.ID %in% training_membership$X1 
  test_membership_file <- paste0("../hapmap3/test-ids-fold", fold, ".txt")
  test_membership <- readr::read_table(test_membership_file, col_names = FALSE)
  test_indic <- gg$fam$family.ID %in% test_membership$X1
  # get column means for snp genotypes matrix of training subjects
  #  tictoc::tic()
  col_means_training <- colMeans(x = gg_imputed[tr_indic, ], na.rm = TRUE)
  col_sds_training <- apply(FUN = sd, X = gg_imputed[tr_indic, ], MARGIN = 2)
  
  GG_centered <- sweep(x = gg_imputed[test_indic, ], 
              MARGIN = 2, 
              STATS = as.array(col_means_training))
  GG <- sweep(x = GG_centered, MARGIN = 2, FUN = "/", STATS = as.array(col_sds_training))
  #tictoc::toc()
  # assign training mean genotype value to any missing values in GG
  #tictoc::tic()
  #tictoc::toc()
  # get indicator that GG SNP is in dbslmm output
  snp_indic <- gg$map$marker.ID %in% dbslmm_output$X1
  
  pgs <- GG[, snp_indic] %*% dbslmm_output$X4
  # store pgs for this fold
  pgs_list[[fold]] <- test_membership %>%
    dplyr::mutate(pgs = pgs[, 1], fold = fold) 
  # extract genotypes for verification set
  verif_indic <- gg$fam$family.ID %in% verif_ids$X1
  # center test genotypes according to mean genotypes in training set
  #tictoc::tic()
  GG_verif_centered <- sweep(x = gg_imputed[verif_indic, ], 
              MARGIN = 2, 
              STATS = as.array(col_means_training))
  GG_verif <- sweep(x = GG_verif_centered, 
                             MARGIN = 2, 
                    FUN = "/",
                             STATS = as.array(col_sds_training))
  muhat_verif <- GG_verif[, snp_indic] %*% dbslmm_output$X4
  muhat[[fold]] <- verif_ids %>%
    dplyr::mutate(muhat = muhat_verif[, 1], fold = fold)
  
} # end loop over folds
# make object phen2 to contain only id columns and a single replicate
phen2 <- test_fam %>%
  dplyr::select(dplyr::all_of((pheno * 5 + 1): (pheno * 5 + 5))) %>%
  as.matrix()
phen <- numeric()
for (row in 1:nrow(phen2)){
  if (isTRUE(all.equal(phen2[row, ] %>% as.numeric(), as.numeric(c(NA, NA, NA, NA, NA))))){
    phen[row] <- NA
  } else {
    foo <- phen2[row, ]
    phen[row] <- foo[!is.na(foo)]
  }
}
# read phen file to get the true phenotype values for the verification set
phen_fn <- paste0("../hapmap3/sim_traits_manual/sims_scenario", scenario, "_", distribution, "_hsq", hsq, ".txt")
ptib <- readr::read_table(phen_fn, col_names = FALSE) %>%
  dplyr::select(1, 2, (pheno + 2)) %>%
  dplyr::right_join(verif_ids, by = c("X1", "X2")) 

#
phen_true_pheno <- test_fam %>%
  dplyr::select(1,2) %>%
  dplyr::mutate(true_pheno = phen) %>%
  dplyr::left_join(ptib, by = c("X1", "X2"))
# get a single vector with true pheno values for all but verification set
true <- numeric()
for (row in 1:nrow(phen_true_pheno)){
  if (isTRUE(all.equal(phen_true_pheno[row, 3:4] %>% as.numeric(), 
                       as.numeric(c(NA, NA))))){
    true[row] <- NA
  } else {
    foo <- phen_true_pheno[row, 3:4]
    true[row] <- foo[!is.na(foo)]
  }
} # end for loop over rows
true_pheno <- phen_true_pheno %>%
  dplyr::select(- c(3,4)) %>%
  dplyr::mutate(true_pheno = true)
pgs_tib <- pgs_list %>%
  # combine rows from pgs_list entries
  dplyr::bind_rows() %>%
  dplyr::arrange(X1) %>%
  dplyr::left_join(true_pheno, by = c("X1", "X2")) %>%
  # calculate residuals
  dplyr::mutate(residual = abs(true_pheno - pgs))
# muhat
n_verif <- nrow(muhat[[1]])
lhs <- numeric()
rhs <- numeric()
n_training <- nrow(pgs_tib)
# handle one verification subject at a time
for (i in 1:n_verif){
  # calculate the arguments to quantile() 
  
  tib2 <- pgs_tib %>%
    dplyr::rename(cv_fold = fold) %>%
    dplyr::mutate(muhat = purrr::map_dbl(.x = cv_fold, 
                                         .f = function(x){
                                           muhat[[x]]$muhat[i]
                                         })) %>%
    dplyr::mutate(muhat_minus_residual = muhat - residual, 
                  muhat_plus_residual = muhat + residual)
  lhs[i] <- quantile(tib2$muhat_minus_residual,
                     probs = floor((n_training + 1) * (alpha / 2)) / n_training)
  rhs[i] <- quantile(tib2$muhat_plus_residual,
                     probs = ceiling((n_training + 1) * (1 - alpha / 2)) / n_training)
}

results <- true_pheno %>%
  dplyr::right_join(verif_ids, by = c("X1", "X2")) %>%
  dplyr::mutate(left_bound = lhs, right_bound = rhs) %>%
  dplyr::mutate(in_interval = true_pheno <= right_bound & true_pheno >= left_bound) %>%
  dplyr::mutate(interval_length = right_bound - left_bound)
interval_lengths <- results %>%
  dplyr::select(interval_length)
coverage <- results %>%  
  dplyr::summarise(coverage = mean(in_interval)) %>%
  unlist()
        
saveRDS(object = coverage, file = paste0("../results/sims-manual-12500-5fold-scenario", scenario, "-", distribution, "-hsq", hsq, "-replicate", pheno, "-alpha", alpha,  ".rds"))
saveRDS(object = interval_lengths, file = paste0("../results/sims-manual-12500-5fold-interval-lengths-scenario", scenario, "-", distribution, "-hsq", hsq, "-replicate", pheno, "-alpha", alpha, ".rds"))

