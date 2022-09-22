#! /usr/bin/env Rscript
rm(list=ls())
library(plyr)
library(bigsnpr)
library(bigreadr)
library(tidyverse)
library(optparse)

# Input parameters
args_list = list(
  make_option("--summ", type="character", default=NULL,
              help="INPUT: gemma file", metavar="character"), 
  make_option("--pathCT", type="character", default=NULL,
              help="INPUT: CT path", metavar="character"),
  make_option("--pathSCT", type="character", default=NULL,
              help="INPUT: SCT  path", metavar="character"),
  make_option("--pheno", type="numeric", default=NULL,
              help="INPUT: phenotype number", metavar="character"), 
  make_option("--cross", type="numeric", default=NULL,
              help="INPUT: cross number", metavar="character"), 
  make_option("--dat", type="character", default=NULL,
              help="INPUT: phenotype type", metavar="character"),
  make_option("--reftype", type="character", default=NULL,
              help="INPUT: reference panel", metavar="character"),
  make_option("--thread", type="numeric", default=1,
              help="OUTPUT: core number", metavar="character"),
  make_option("--scenario", type = "character", default = NULL, metavar = "character"),
  make_option("--distribution", type = "character", default = NULL, metavar = "character"),
  make_option("--hsq", type = "numeric", default = NULL, metavar = "character")
  
  
)

opt_parser = OptionParser(option_list=args_list)
opt = parse_args(opt_parser)

# opt <- list(summ = "/net/mulan/disk2/yasheng/comparisonProject/05_internal_c/pheno1/output/summary_ukb_cross1.assoc.txt",
#             path = "/net/mulan/disk2/yasheng/comparisonProject/05_internal_c/pheno1/",
#             thread = 4,
#             pheno = 1,
#             cross = 1,
#             dat = "continuous",
#             reftype = "ukb")

# set parameter
#comp_str <- "/net/mulan/disk2/yasheng/comparisonProject/"
if (opt$dat == "binary"){
  ref_str <- paste0(comp_str, "04_reference/", opt$reftype, 
                  "/geno_inter/merge_b")
  val_str <- paste0(comp_str, "03_subsample/", opt$dat, "/pheno", opt$pheno,
                    "/val/impute_inter/merge")
}
if (opt$dat == "continuous"){
#  ref_str <- paste0(comp_str, "04_reference/", opt$reftype,
#                    "/geno_inter/merge_c")
  ref_str <- paste0("../dat-quant/reference/chr1")
#  val_str <- paste0(comp_str, "03_subsample/", opt$dat, "/pheno", opt$pheno,
#                    "/val/",opt$reftype,"/impute_inter/merge")
  val_str <- paste0("../dat-quant/validation/chr1")
}

ref_sub_str <- paste0(ref_str, "_sub-", as.numeric(as.POSIXlt(Sys.time())))
val_sub_str <- paste0(val_str, "_sub-", as.numeric(as.POSIXlt(Sys.time())))





# interect val and ref data
# start <- proc.time()
ref_bed <- paste0(ref_str, ".rds")
val_bed <- paste0(val_str, ".rds")



#ref_bed <- snp_attach(paste0(ref_str, ".rds"))
ref_bed <- snp_attach(ref_bed)
val_bed <- snp_attach(val_bed)
cat("ref and val are loaded!\n")
# val_snp <- fread2(paste0(val_str, ".bim"))
val_bed$genotypes <- snp_fastImputeSimple(val_bed$genotypes) #impute missing values
# process ref data
ref_G <- snp_fastImputeSimple(ref_bed$genotypes) # impute missing values
ref_CHR <- ref_bed$map$chromosome
ref_POS <- ref_bed$map$physical.pos
ref_n_snp <- dim(ref_G)[2]
ref_map <- ref_bed$map[, -3]
names(ref_map) <- c("chr", "rsid", "pos", "a1", "a0")


# process summary statistics
summstats <- fread2(opt$summ, select =  c(1, 2, 3, 7, 6, 9, 10))
colnames(summstats) <- c("chr", "rsid", "pos", "a0", "a1", "beta", "se") # calculate P
t <- summstats$beta/summstats$se
p_val <- ifelse(t < 0, pnorm(t), pnorm(t, lower.tail = F))*2
summstats$pval <- ifelse(p_val == 0,
                         min(p_val[-which(p_val==0)]),
                         p_val)

# map summary statistics and reference panel
info_snp <- snp_match(summstats, ref_map)
beta <- rep(0, ref_n_snp)
lp_val <- rep(0, ref_n_snp)
beta[ref_map[, 2]%in%info_snp[, 5]] <- info_snp$beta
lp_val[ref_map[, 2]%in%info_snp[, 5]] <- -log10(info_snp$pval)

# clump
all_keep <- snp_grid_clumping(ref_G,
                              # grid.thr.r2 = c(0.8),
                              # grid.base.size = 500,
                              infos.chr = ref_CHR,
                              infos.pos = ref_POS,
                              lpS = lp_val,
                              ncores = opt$thread)
cat("Clumping is finished!\n")



# threshold
ct_bk_str <- ifelse(opt$dat == "continuous",
                    paste0("../dat-quant/CT/scenario", opt$scenario, "/", opt$distribution, "/hsq", opt$hsq, 
                           "/summary_cross", opt$cross, "_sub-",
                           as.numeric(as.POSIXlt(Sys.time()))),
                    paste0(comp_str, "06_internal_b/pheno", opt$pheno, 
                           "/CT/summary_cross", opt$cross, "_sub-",
                           as.numeric(as.POSIXlt(Sys.time()))))

val2_bed <- val_bed %>%
  snp_subset(ind.col = rows_along(ref_bed$map))
val2 <- snp_attach(val2_bed)


multi_PRS <- snp_grid_PRS(val2$genotypes,
                          all_keep,
                          betas = beta,
                          lpS = lp_val,
                          # n_thr_lpS = 2,
                          backingfile = ct_bk_str,
                          ncores = opt$thread)
nn <- nrow(attr(all_keep, "grid"))
grid2 <- attr(all_keep, "grid") %>%
  mutate(thr.lp = list(attr(multi_PRS, "grid.lpS.thr")), id = c(1:nn)) %>%
  unnest(cols = "thr.lp")
s <- nrow(grid2)
cat("Threashold is finished!\n")

# load validation trait
if (opt$dat == "continuous"){
#  y <- fread2(paste0(comp_str, "03_subsample/continuous/pheno", opt$pheno, 
#                     "/val/",opt$reftype,"/02_pheno_c.txt"))[, 1]
  y <- fread2(paste0("../dat-quant/validation/scenario", opt$scenario, "_", 
                    opt$distribution, "_hsq", opt$hsq, "_replicate", opt$pheno, ".txt"))[, 1]
} else {
  y <- fread2(paste0(comp_str, "03_subsample/binary/pheno", opt$pheno, 
                     "/val/02_pheno_b.txt"))[, 1]
  covar <- fread2(paste0(comp_str, "03_subsample/binary/pheno", opt$pheno, 
                         "/val/03_cov_eff.txt"))[, 1]
}

val_bk_str2 <- paste0(val_str, "_sub-", as.numeric(as.POSIXlt(Sys.time())))
if (all(is.na(y)) == T){
  idx <- which(!is.na(y))
  val_bed <- snp_attach(snp_subset(val_bed, ind.row = idx, 
                                   backingfile = val_bk_str2))
  y <- y[idx]
} 

## subsample phenotype
if (opt$dat == "continuous"){
  
  grid2$valIdx <- big_apply(multi_PRS, a.FUN = function(X, ind, s, y.train) {
    #single_PRS <- rowSums(X[, ind + s * (0:21)])
    single_PRS <- X[, ind]
    return(cor(single_PRS, y.train)^2)
  },
  ind = 1:s,
  s = s,
  y.train = y,
  a.combine = 'c',
  block.size = 1,
  ncores = opt$thread
  )
  all_pgs_CT <- big_apply(multi_PRS, a.FUN = function(X, ind, s, y.train) {
    #single_PRS <- rowSums(X[, ind + s * (0:21)])
    single_PRS <- X[, ind]
    return(single_PRS)
  },
  ind = 1:s,
  s = s,
  y.train = y,
  a.combine = 'cbind',
  block.size = 1,
  ncores = opt$thread
  )
  final_mod <- snp_grid_stacking(multi_PRS,
                                 y,
                                 ncores = opt$thread,
                                 K = 10)
} else{
  
  grid2$valIdx <- big_apply(multi_PRS, a.FUN = function(X, ind, s, y.train, covar.train) {
    single_PRS <- rowSums(X[, ind + s * (0:21)]+covar.train)
    return(bigstatsr::AUC(single_PRS, y.train))
  },
  ind = 1:s,
  s = s,
  y.train = y,
  covar.train = as.vector(covar),
  a.combine = 'c',
  block.size = 1,
  ncores = opt$thread
  )
  final_mod <- snp_grid_stacking(multi_PRS,
                                 y,
                                 covar.train = as.matrix(covar),
                                 ncores = opt$thread,
                                 K = 10)
}

# snp_info_CT <- data.frame(ref_map$rsid, ref_map$a1)

# idx_mat_CT <- t(plyr::aaply(c(1: 1400), 1, function(ss){
#   info_prs <- grid2 %>% slice(ss)
#   c_idx <- c(1: n_snp) %in% unlist(map(all_keep, info_prs$id))
#   t_idx <- c(1: n_snp) %in% which(lp_val >= 0.999999*info_prs$thr.lp)
#   return(ifelse(c_idx==T&t_idx==T, T, F))
# }))

# save(snp_info_CT, idx_mat_CT, all_pgs_CT,
#      file = paste0(opt$path, "CT/res_",
#                    opt$reftype, "_cross", opt$cross, ".RData"))

# CT
max_prs <- grid2 %>% arrange(desc(valIdx)) %>% slice(1)
c_idx <- c(1: ref_n_snp) %in% unlist(map(all_keep, max_prs$id))
t_idx <- c(1: ref_n_snp) %in% which(lp_val >= 0.999999*max_prs$thr.lp)
idx <- ifelse(c_idx==T&t_idx==T, T, F)
snp_sig_CT <- data.frame(ref_map$rsid[idx],
                         ref_map$a1[idx],
                         beta[idx])

# SCT
new_beta <- final_mod$beta.G 
idx <- which(new_beta != 0)
snp_sig_SCT <- data.frame(ref_map$rsid[idx],
                          ref_map$a1[idx],
                          new_beta[idx])

# # output
write.table(grid2, file = paste0(opt$pathCT, "grid2_", opt$reftype,  "_pheno", opt$pheno, "_cross", opt$cross, ".txt"),
             col.names = F, row.names = F, quote = F)
write.table(snp_sig_CT, file = paste0(opt$pathCT, "esteff_", opt$reftype,  "_pheno", opt$pheno, "_cross", opt$cross, ".txt"),
             col.names = F, row.names = F, quote = F)
write.table(snp_sig_SCT, file = paste0(opt$pathSCT, "esteff_", opt$reftype, "_pheno", opt$pheno, 
                                       "_cross", opt$cross, ".txt"),
             col.names = F, row.names = F, quote = F)

# remove file
system(paste0("rm ", ct_bk_str, ".bk"))
system(paste0("rm ", ct_bk_str, ".rds"))
# system(paste0("rm ", ref_sub_str, ".bk"))
# system(paste0("rm ", ref_sub_str, ".rds"))
# system(paste0("rm ", val_sub_str, ".bk"))
# system(paste0("rm ", val_sub_str, ".rds"))
if(file.exists(paste0(val_bk_str2, ".bk"))){
  system(paste0("rm ", val_bk_str2, ".bk"))
  system(paste0("rm ", val_bk_str2, ".rds"))
}
