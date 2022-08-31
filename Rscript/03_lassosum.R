#! /usr/bin/env Rscript
rm(list = ls())
library(plyr)
library(bigreadr)
library(optparse)
library(lassosum)
library(doParallel)

# Parameter setting
args_list = list(
  make_option("--summ", type="character", default=NULL,
              help="INPUT: summary data prefix", metavar="character"), 
  make_option("--pheno", type="numeric", default=NULL,
              help="INPUT: phenotype number", metavar="character"),
  make_option("--dat", type="character", default="ukb",
              help="INPUT: input reference panel type", metavar="character"),
  make_option("--n", type="numeric", default=NULL,
              help="INPUT: sample size", metavar="character"),
  make_option("--thread", type="numeric", default=1,
              help="INPUT: thread", metavar="character"),
  make_option("--cross", type="numeric", default=NULL,
              help="INPUT: cross", metavar="character"),
  make_option("--reftype", type="character", default="ukb",
              help="INPUT: input reference panel type", metavar="character"),
  make_option("--path", type="character", default=NULL,
              help="OUTPUT: output path", metavar="character")
)

opt_parser = OptionParser(option_list=args_list)
opt = parse_args(opt_parser) 

# opt <- list(summ = "/net/mulan/disk2/yasheng/comparisonProject/05_internal_c/pheno1/output/summary_ukb_cross1.assoc.txt",
#             thread = 4,
#             pheno = 1,
#             n = 268505,
#             cross = 1,
#             dat = "continuous",
#             reftype = "ukb",
#             path = "/net/mulan/disk2/yasheng/comparisonProject/05_internal_c/pheno1/lassosum")

# parameters
comp_str <- "/net/mulan/disk2/yasheng/comparisonProject/"

if (opt$dat == "continuous"){
  ref_str <- paste0(comp_str, "04_reference/", opt$reftype,
                    "/geno_inter/merge_c")
  valid_str <- paste0(comp_str, "03_subsample/", opt$dat, "/pheno",
                      opt$pheno,"/val/",opt$reftype,"/impute_inter/merge")
  
} else {
  ref_str <- paste0(comp_str, "04_reference/", opt$reftype,
                    "/geno_inter/merge_b")
  valid_str <- paste0(comp_str, "03_subsample/", opt$dat, "/pheno",
                      opt$pheno,"/val/impute_inter/merge")
  
}

# summary statistics input
summstats <- fread2(opt$summ, select =  c(1, 2, 3, 6, 7, 9, 10))
colnames(summstats) <- c("chr", "snp", "pos", "a1", "a2", "beta", "se") # calculate P 
summstats <- summstats[!is.nan(summstats$se), ]

# p to cor
t <- summstats$beta/summstats$se
p_val <- ifelse(t < 0, pnorm(t), pnorm(t, lower.tail = F))*2
p_val_z <- ifelse(p_val == 0, min(p_val[-which(p_val==0)]), p_val)
cor <- p2cor(p = p_val_z, n = opt$n, sign = summstats$beta)

# block information
setwd(system.file("data", package="lassosum"))
LDblocks <- "EUR.hg19"

# parallel setting
threads <- as.numeric(opt$thread)
cl <- makeCluster(threads, type="FORK")

# estimation
out <- lassosum.pipeline(cor = cor, chr = summstats$chr, pos = summstats$pos,
                         A1 = summstats$a1, A2 = summstats$a2, destandardize = F,
                         ref.bfile = ref_str, test.bfile = ref_str,
                         LDblocks = LDblocks, cluster = cl)

# validation
if (opt$dat == "continuous"){
  idx_val <- fread2(paste0(comp_str, "03_subsample/continuous/pheno",opt$pheno,"/val/",opt$reftype,"/01_idx.txt"))
  pheno_val_tot <- fread2(paste0(comp_str, "03_subsample/continuous/pheno",
                                 opt$pheno, "/val/",opt$reftype,"/02_pheno_c.txt"))[, 1]
  pheno_val <- data.frame(FID = idx_val[, 1], IID = idx_val[, 1],
                          pheno_val_tot)
  valid <- validate(out, pheno = pheno_val, plot = F, 
                    test.bfile = valid_str)
}
if (opt$dat == "binary"){
  cat("binary cov\n")
  idx_val <- fread2(paste0(comp_str, "03_subsample/binary/pheno",opt$pheno,"/val/01_idx.txt"))
  pheno_val_tot <- fread2(paste0(comp_str, "03_subsample/binary/pheno",
                                 opt$pheno, "/val/02_pheno_b.txt"))[, 1]
  pheno_val <- data.frame(FID = idx_val[, 1], IID = idx_val[, 1],
                          pheno_val_tot)[!is.na(pheno_val_tot), ]
  covar <- fread2(paste0(comp_str, "03_subsample/binary/pheno", opt$pheno, 
                         "/val/03_cov_eff.txt"))[, 1]
  cov_dat <- data.frame(FID = idx_val[, 1], IID = idx_val[, 1], 
                        covar)
  valid <- validate(out, pheno = pheno_val, covar = cov_dat, plot = F, 
                    test.bfile = valid_str)
}



# snp_info_lassosum <- data.frame(snp = summstats$snp[out$sumstats$order],
#                                 a1 = summstats$a1[out$sumstats$order])
# beta_lassosum <- out$beta
# pgs_lassosum <- out$pgs
# save(snp_info_lassosum, beta_lassosum, pgs_lassosum,
#      file = paste0(opt$path, "/res_", opt$reftype, "_cross", opt$cross, ".RData"))

# write.table(valid$validation.table,
#             file = paste0(opt$path, "/val_", opt$reftype,
#                           "_cross", opt$cross, ".txt"),
#             row.names = F, col.names = F, quote = F)

# best
out_best <- subset(out, s = valid$best.s, lambda = valid$best.lambda)
esteff_best <- data.frame(snp = summstats$snp[out_best$sumstats$order],
                          a1 = summstats$a1[out_best$sumstats$order],
                          beta = out_best$beta)
esteff_best_nz <- esteff_best[esteff_best[, 3] != 0, ]

# write.table(esteff_best_nz,
#             file = paste0(opt$path, "/esteff_", opt$reftype,
#                           "_cross", opt$cross, ".txt"),
#             row.names = F, col.names = F, quote = F)
