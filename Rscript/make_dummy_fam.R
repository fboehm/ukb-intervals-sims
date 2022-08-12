library(magrittr)

fam_file <- snakemake@input[[1]]
vroom::vroom(fam_file, col_names = FALSE) %>%
  dplyr::select(1:6) %>%
  vroom::vroom_write(file = snakemake@output[[1]], col_names = FALSE)
