# Goal is to read in the HapMap3 plink map file and to output a plain text file with all
# rsids for the HM3 SNPs

library(magrittr)
map <- vroom::vroom("../hapmap3/hapmap3_r3_b36_fwd.consensus.qc.poly.map", col_names = FALSE)
map %>%
  dplyr::select(X2) %>% # choose the one column with rsids
  vroom::vroom_write("../hapmap3/snp_list.txt", col_names = FALSE)
