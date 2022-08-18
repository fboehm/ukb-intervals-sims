library(magrittr)

#fam_file <- snakemake@input[[1]]
fam_file <- "/net/mulan/disk2/yasheng/predictionProject/plink_file/ukb/chr22.fam"
vroom::vroom(fam_file, col_names = FALSE) %>%
  dplyr::select(1:4) %>%
  dplyr::mutate(X5 = NA, X6 = NA) %>%
  #vroom::vroom_write(file = snakemake@output[[1]], col_names = FALSE)
  vroom::vroom_write(file = "../dat/chr22.fam", col_names = FALSE)
