# Goal: use DBSLMM-estimated SNP effects and genotypes to calculate PGS
args = commandArgs(trailingOnly=TRUE)
(scenario <- args[1])
(distribution <- args[2])
(hsq <- as.numeric(args[3]))
(pheno <- as.numeric(args[4]))
library(magrittr)
#scenario_vec <- c("I", "II", "III", "IV")
#distribution_vec <- c("laplace", "normal", "scaledt")
#hsq_vec <- c(0.1, 0.2, 0.5)
#settings_tib <- expand.grid(hsq_vec, distribution_vec, scenario_vec) %>%
#  tibble::as_tibble() %>%
#  dplyr::rename(hsq = 1, distribution = 2, scenario = 3)


#for (scenario in scenario_vec){
#  for (distribution in distribution_vec){
#    for (hsq in hsq_vec){



#myfun <- function(scenario, distribution, hsq){
      n_folds <- 5
      alpha <- 0.1
      b_file <- "../hapmap3/plink_files_for_sims/chr1.rds"
      gg <- bigsnpr::snp_attach(b_file)
      # read ids for verification set
      v_file <- "../hapmap3/subjects_for_sims_verification.txt"
      verif_ids <- readr::read_table(v_file, col_names = FALSE)

      # read phen file with "true" phenotypes
      phen_file <- paste0("../hapmap3/sim_traits/sims_scenario", scenario, "_", distribution, "_hsq", hsq, ".txt.phen")
      phen <- readr::read_table(phen_file, col_names = FALSE)
      coverage <- numeric() 
#      tictoc::tic()
      #for (pheno in 1:10){
        pgs_list <- list()
        muhat <- list() # to store the muhat (Xn+1) values from Eqn 3.1
        for (fold in 1:n_folds){
          # DBSLMM outputted file
          d_file <- paste0("../dat-quant/DBSLMM/scenario", scenario, "/", distribution, "/hsq", hsq, "/summary_ukb_pheno", pheno, "_scenario", scenario, "_", distribution, "_hsq", hsq, "_fold", fold, "_chr1_best.dbslmm.txt")
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
          col_means_training <- colMeans(x = gg$genotypes[tr_indic, ], na.rm = TRUE)
         # tictoc::toc()
        #  tictoc::tic()
#          col_means_trainingA <- bigstatsr::big_apply(X = gg$genotypes, 
#                                                     a.FUN = function(X, ind, row_indicator){
#                                                       colMeans(X[row_indicator, ind], na.rm = TRUE)
#                                                     }, row_indicator = tr_indic, 
#                                                     a.combine = 'c', 
#                                                     ncores = bigstatsr::nb_cores())
        #  tictoc::toc()
          # center test genotypes according to mean genotypes in training set
         # tictoc::tic()
#          GGA <- bigstatsr::big_apply(X = gg$genotypes, 
#                                     a.FUN = function(X, ind, row_indicator, col_means_training){
#                                       sweep(x = X[row_indicator, ind], MARGIN = 2, STATS = as.array(col_means_training[ind]))
#                                     },
#                                     a.combine = 'cbind',
#                                     row_indicator = test_indic, col_means_training = col_means_trainingA,
#                                     ncores = bigstatsr::nb_cores())
          #tictoc::toc()
          #tictoc::tic()
          GG <- sweep(x = gg$genotypes[test_indic, ], 
                      MARGIN = 2, 
                      STATS = as.array(col_means_training))
          #tictoc::toc()
          # assign training mean genotype value to any missing values in GG
          #tictoc::tic()
          for (col in 1:ncol(GG)){
            foo <- GG[, col]
            foo[is.na(foo)] <- 0 #col_means_training[col]
            GG[, col] <- foo
          }
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
          GG_verif <- sweep(x = gg$genotypes[verif_indic, ], 
                      MARGIN = 2, 
                      STATS = as.array(col_means_training))
          #tictoc::toc()
          # assign training mean genotype value to any missing values in GG_verif
          #tictoc::tic()
          for (col in 1:ncol(GG_verif)){
            foo <- GG_verif[, col]
            foo[is.na(foo)] <- 0 #col_means_training[col]
            GG_verif[, col] <- foo
          }
          #tictoc::toc()
          muhat_verif <- GG_verif[, snp_indic] %*% dbslmm_output$X4
          muhat[[fold]] <- verif_ids %>%
            dplyr::mutate(muhat = muhat_verif[, 1], fold = fold)
          
        } # end loop over folds
                # make object phen2 to contain only id columns and a single replicate
        phen2 <- phen %>%
          dplyr::select(1, 2, dplyr::all_of(2 + pheno)) %>%
          dplyr::rename(true_pheno = 3)
        pgs_tib <- pgs_list %>%
          # combine rows from pgs_list entries
          dplyr::bind_rows() %>%
          dplyr::arrange(X1) %>%
          dplyr::left_join(phen2, by = c("X1", "X2")) %>%
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

        coverage <- phen2 %>%
          dplyr::right_join(verif_ids, by = c("X1", "X2")) %>%
          dplyr::mutate(left_bound = lhs, right_bound = rhs) %>%
          dplyr::mutate(in_interval = true_pheno <= right_bound & true_pheno >= left_bound) %>%
          dplyr::summarise(coverage = mean(in_interval)) %>%
          unlist()
        
        
        
#      } # end loop over replicates
#      tictoc::toc()
#      return(coverage)
#}      
#    }
#  }
#}
#library(future)
#plan(multicore, workers = 36)
#tictoc::tic()

#coverages <- settings_tib %>%
#  dplyr::mutate(coverage = furrr::future_pmap(.l = settings_tib, .f = myfun))
#tictoc::toc()
saveRDS(object = coverage, file = paste0("../results/sims-12500-5fold-scenario", scenario, "-", distribution, "-hsq", hsq, "-replicate", pheno, ".rds"))
