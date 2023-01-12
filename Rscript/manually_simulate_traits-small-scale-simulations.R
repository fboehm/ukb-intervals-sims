# GOAL: Manually simulate traits from SNP effects files

library(magrittr)

# Read genotypes from bed file and center & standardize each SNP within the training set
bed_file <- "../hapmap3-small-scale-simulations/plink_files_for_sims/chr1.bed"
bigsnpr::snp_readBed(bed_file)
rds_file <- "../hapmap3-small-scale-simulations/plink_files_for_sims/chr1.rds"
gg <- bigsnpr::snp_attach(rds_file)
gg2 <- bigsnpr::snp_fastImputeSimple(Gna = gg$genotypes)
gg_scaled <- bigstatsr::big_apply(X = gg2, 
                                  a.FUN = function(X, ind){
                                      apply(FUN = scale, X = X[, ind], MARGIN = 2)
                                  },
                                  a.combine = "cbind")

scenario_vec <- c("I", "II", "III", "IV")
distribution_vec <- c("laplace", "normal", "scaledt")
hsq_vec <- c(0.1, 0.2, 0.5)
## SCENARIO I
# normal distribution
output_dir <- "../hapmap3-small-scale-simulations/sim_traits_manual/"
if (!dir.exists(output_dir)){
    dir.create(output_dir, recursive = TRUE)
}
set.seed(2023-01-11)
for (scenario in scenario_vec){
    for (distribution in distribution_vec){
        for (hsq in hsq_vec){
            infile <- paste0("../hapmap3/snp_effects/scenario", scenario, "_", distribution, "_hsq", hsq, ".txt")
            snp_effects <- vroom::vroom(file = infile, col_names = FALSE)
            xb_pre <- gg_scaled %*% snp_effects$X2
            xb <- scale(xb_pre) * sqrt(hsq)
            y <- matrix(nrow = length(xb), ncol = 10)
            for (replicate in 1:10){
                eps_pre <- rnorm(n = length(xb), mean = 0, sd = 1)
                eps <- scale(eps_pre) * sqrt(1 - hsq)
                y[, replicate] <- scale(xb + eps)
            }
            
            cbind(gg$fam[, 1:2], y) %>%
                tibble::as_tibble() %>%
                vroom::vroom_write(file = paste0(output_dir, "sims_scenario", scenario, "_", distribution, "_hsq", hsq, ".txt"), col_names = FALSE)
        }
    }
}
