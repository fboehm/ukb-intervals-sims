#! /usr/bin/env Rscript
rm(list=ls())
library(bigstatsr)
library(plyr)
library(bigsnpr)
library(bigreadr)
library(optparse)
library(tidyverse)

## Input parameters
args_list = list(
  make_option("--summ", type="character", default=NULL,
              help="INPUT: gemma file", metavar="character"), 
  make_option("--LDpred2Path", type="character", default=NULL,
              help="INPUT: CT path", metavar="character"),
  make_option("--pheno", type="numeric", default=NULL,
              help="INPUT: phenotype number", metavar="character"), 
  make_option("--cross", type="numeric", default=NULL,
              help="INPUT: cross number", metavar="character"), 
  make_option("--dat", type="character", default=NULL,
              help="INPUT: phenotype type", metavar="character"),
  make_option("--reftype", type="character", default=NULL,
              help="INPUT: phenotype type", metavar="character"),
  make_option("--thread", type="numeric", default=1,
              help="OUTPUT: core number", metavar="character"), 
  make_option("--chr", type="numeric", default=NULL,
              help="OUTPUT: chr number", metavar="character")
  ,
  make_option("--model", type="character", default=NULL,
              help="INPUT: model setting", metavar="character")
)

opt_parser = OptionParser(option_list=args_list)
opt = parse_args(opt_parser)

# opt <- list(summ = "/net/mulan/disk2/yasheng/comparisonProject/06_internal_b/pheno1/output/summary_rare_cross1_chr",
#             LDpred2Path = "/net/mulan/disk2/yasheng/comparisonProject/06_internal_b/pheno15/LDpred2/",
#             thread = 1,
#             pheno = 1,
#             cross = 1,
#             dat = "binary",
#             reftype = "hm3_rare",
#             chr = 22)


# set parameter
p_len <- 21 # length of vector of p's
#comp_str1 <- "/net/mulan/disk2/yasheng/comparisonProject/"
#comp_str2 <- "/net/mulan/disk2/yasheng/comparisonProject/"
comp_str1 <- "~/research/ukb-intervals-sims/dat/"
comp_str2 <- comp_str1

if (opt$dat == "binary"){
  ref_str <- paste0(comp_str2, "04_reference/", opt$reftype, 
                  "/geno_inter/chr", opt$chr, "_b")
  val_str <- paste0(comp_str1, "03_subsample/", opt$dat, "/pheno", opt$pheno, 
                    "/val/impute_inter/chr", opt$chr)
} else {
  #ref_str <- paste0(comp_str2, "04_reference/", opt$reftype, 
  #                  "/geno_inter/chr", opt$chr, "_c")
  # define path to chr-specific plink files for ref data
  ref_str <- paste0(comp_str1, "reference/ukb/geno/chr", opt$chr)
  # define path to chr-specific validation plink files
  val_str <- paste0(comp_str1, "validation/geno/chr", opt$chr)
  #val_str <- paste0(comp_str1, "03_subsample/", opt$dat, "/pheno", opt$pheno, 
  #                  "/val/", opt$reftype, "/impute_inter/chr", opt$chr)
}
ref_sub_str <- paste0(ref_str, "_sub-", as.numeric(as.POSIXlt(Sys.time())))
# val_str <- paste0(comp_str1, "03_subsample/", opt$dat, "/pheno", opt$pheno, 
#                   "/val/rare/impute_inter/chr", opt$chr)
val_sub_str <- paste0(val_str, "_sub-", as.numeric(as.POSIXlt(Sys.time())))

# interect val and ref data
if(!file.exists(paste0(val_str, ".rds")) | !file.exists(paste0(val_str, ".bk"))){
  # system(paste0("rm ", val_str, ".bk"))
  if(file.exists(paste0(val_str, ".bk"))){
    system(paste0("rm ", val_str, ".bk"))
  }
  val_bed <- snp_readBed(paste0(val_str, ".bed"))
}
if(!file.exists(paste0(ref_str, ".rds")) | !file.exists(paste0(ref_str, ".bk"))){
  # system(paste0("rm ", val_str, ".bk"))
  if(file.exists(paste0(ref_str, ".bk"))){
    system(paste0("rm ", ref_str, ".bk"))
  }
  ref_bed <- snp_readBed(paste0(ref_str, ".bed"))
}

ref_bed <- snp_attach(paste0(ref_str, ".rds"))
val_bed <- snp_attach(paste0(val_str, ".rds"))

length(ref_bed$map$marker.ID)
length(val_bed$map$marker.ID)

if(all(ref_bed$map$marker.ID == val_bed$map$marker.ID) == F){

  snp_inter <- intersect(val_bed$map$marker.ID, ref_bed$map$marker.ID)
  ref_bed <- snp_attach(snp_subset(ref_bed, 
                                   ind.col = which(ref_bed$map$marker.ID%in%snp_inter), 
                                   backingfile = ref_sub_str))
  val_bed <- snp_attach(snp_subset(val_bed, 
                                   ind.col = which(val_bed$map$marker.ID%in%snp_inter), 
                                   backingfile = val_sub_str))
}

# process ref data
# ref_bed <- snp_readBed(paste0(ref_str, ".bed"))
ref_bed$map$genetic.dist <- snp_asGeneticPos(infos.chr = ref_bed$map$chromosome, 
                                             infos.pos = ref_bed$map$physical.pos)
ref_G <- ref_bed$genotypes
ref_CHR <- ref_bed$map$chromosome
ref_POS <- ref_bed$map$physical.pos
ref_n_snp <- dim(ref_G)[2]
ref_map <- ref_bed$map
names(ref_map) <- c("chr", "rsid", "g.dis", "pos", "a1", "a0")

## process summary data
summstats <- fread2(paste0(as.character(opt$summ), opt$chr, ".assoc.txt"),
                    select =  c(1, 2, 3, 5, 6, 7, 9, 10))
colnames(summstats) <- c("chr", "rsid", "pos", "n_eff", 
                         "a1", "a0",  "beta", "beta_se")
summstats <- summstats[!is.nan(summstats$beta_se),]
info_snp <- snp_match(summstats, ref_map[, -c(2, 3)], match.min.prop = 0.05)
df_beta <- info_snp[, c("beta", "beta_se", "n_eff")]
info_chr <- info_snp$`_NUM_ID_`
POS2 <- ref_map$g.dis[ref_map$rsid%in%info_snp$rsid]

# est h2
#if (opt$reftype == "hm3"){
if (opt$reftype == "ukb"){
  corr <- snp_cor(ref_G, ind.col = info_chr, 
                  ncores = opt$thread,
                  infos.pos = POS2,
                  size = 3 / 1000)
} else {
  corr <- snp_cor(ref_G, ind.col = info_chr, 
                  ncores = opt$thread,
                  size = 200)
}

ldsc <- snp_ldsc2(corr, df_beta)
h2_est <- ldsc[["h2"]]
cat (opt$chr, ":", h2_est, "\n")

if(h2_est < 0){
  beta_LDpred2 <- NA
} else {
  corr_sp <- bigsparser::as_SFBM(as(corr, "dgCMatrix"))
  if (opt$model == "LDpred2-inf"){
  beta_inf <- snp_ldpred2_inf(corr_sp, df_beta, h2 = h2_est)
  cat ("Inf model is ok!\n")
  }
  
  ## LDpred2-auto
  if (opt$model == "LDpred2-auto"){
  auto <- snp_ldpred2_auto(corr_sp, df_beta, h2_init = h2_est)
  beta_auto <- auto[[1]]$beta_est
  beta_auto <- ifelse(is.na(beta_auto), 0, beta_auto)
  cat ("LDpred-auto model is ok!\n")
  }
  
  ## LDpred2
  if (opt$model == "LDpred2-m"){
  
  
  if (opt$dat == "continuous"){
    #y <- fread2(paste0(comp_str1, "03_subsample/continuous/pheno", opt$pheno, 
    #                   "/val/", opt$reftype, "/02_pheno_c.txt"))[, 1]
    y <- vroom::vroom(paste0(comp_str1, "hsq", hsq, "_pcausal", pcausal, "/validation/pheno", opt$pheno, "_hsq", hsq, "_pcausal", pcausal, ".txt"),
                      col_names = FALSE)[, 1]
  } else {
    y <- fread2(paste0(comp_str1, "03_subsample/binary/pheno", opt$pheno, 
                       "/val/02_pheno_b.txt"))[, 1]
    covar <- fread2(paste0(comp_str1, "03_subsample/binary/pheno", opt$pheno, 
                           "/val/03_cov_eff.txt"))[, 1]
  }
  
  # p_seq <- signif(seq_log(1e-5, 1, length.out = 2), 2)
  # h_seq <- round(h2_est * 1, 4)
  
  p_seq <- signif(seq_log(1e-5, 1, length.out = p_len), 2)
  h_seq <- round(h2_est * c(0.7, 1, 1.4), 4)
  params <- expand.grid(p = p_seq, h2 = h_seq, sparse = c(FALSE, TRUE))
  beta_grid <- snp_ldpred2_grid(corr_sp, df_beta,
                                params, ncores = opt$thread)
  val_sub_str2 <- paste0(val_str, "_sub-", as.numeric(as.POSIXlt(Sys.time())))
  if (any(is.na(y)) == T){
    idx <- which(!is.na(y))
    val_bed <- snp_attach(snp_subset(val_bed, 
                                     ind.row = idx, 
                                     backingfile = val_sub_str2))
    y <- y[idx]
  }
  
  val_G <- val_bed$genotypes
  pred_grid <- big_prodMat(val_G, beta_grid, ind.col = info_chr)
  idx_na <- apply(pred_grid, 2, function(a) all(is.na(a))) | 
    apply(beta_grid, 2, function(a) all(is.na(a)))
  beta_grid_na <- beta_grid[, !idx_na]
  pred_grid_na <- pred_grid[, !idx_na]
  params_na <- params[!idx_na, ]
   save(beta_grid, pred_grid, file = paste0(opt$LDpred2Path, "res_", 
                                            opt$ref, "_pheno", opt$pheno, 
                                            "_cross", opt$cross, 
                                            "_chr", opt$chr, ".RData"))
  
  if (all(idx_na)){
    
    beta_LDpred2 <- data.frame(info_snp$rsid,
                               info_snp$a1,
                               beta_inf,
                               0,
                               0,
                               beta_auto)
  } else {
    
    
    if (opt$dat == "continuous"){
      params_na[c("coef", "score")] <-
        big_univLinReg(big_copy(pred_grid_na), y)[c("estim", "score")]
      params_na$idx_val <- apply(pred_grid_na, 2, function(a) cor(a, y)^2)
    } else {
      params_na[c("coef", "score")] <-
        big_univLinReg(big_copy(pred_grid_na), 
                       y, 
                       covar.train = as.matrix(covar))[c("estim", "score")]
      params_na$idx_val  <- apply(pred_grid_na, 2, AUC, target = y)
    }
    
    ## parameters
    params_na %>%
      mutate(sparsity = colMeans(beta_grid_na == 0), id = c(1: nrow(params_na))) %>%
      arrange(desc(score)) %>%
      mutate_at(4:8, signif, digits = 3)
    # no-sparsity effect
    best_grid_nosp <- params_na %>%
      mutate(id = c(1: nrow(params_na))) %>%
      filter(!sparse) %>%
      arrange(desc(score)) %>%
      slice(1) %>%
      { beta_grid_na[, .$id] * .$coef }
    # sparsity effect
    best_grid_sp <- params_na %>%
      mutate(id = c(1: nrow(params_na))) %>%
      filter(sparse) %>%
      arrange(desc(score)) %>%
      slice(1) %>%
      { beta_grid_na[, .$id] * .$coef }
    cat ("LDpred model is ok!\n")
    ## output
    beta_LDpred2 <- data.frame(info_snp$rsid,
                                info_snp$a1,
                                beta_inf,
                                best_grid_nosp,
                                best_grid_sp,
                                beta_auto)
    }
  }
}

 write.table(beta_LDpred2,
             file = paste0(opt$LDpred2Path, "esteff_", opt$reftype, "_pheno", opt$pheno,
                           "_cross", opt$cross, "_chr", opt$chr, ".txt"),
             col.names = F, row.names = F, quote = F)

system(paste0("rm ", val_sub_str2, ".bk"))
system(paste0("rm ", val_sub_str, ".bk"))
system(paste0("rm ", val_sub_str, ".rds"))
system(paste0("rm ", ref_sub_str, ".bk"))
system(paste0("rm ", ref_sub_str, ".rds"))

