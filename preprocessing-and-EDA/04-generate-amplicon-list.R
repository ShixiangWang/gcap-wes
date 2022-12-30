setwd(file.path(PROJ_DIR, "preprocessing-and-EDA/"))

library(magrittr)

load("data/pancan_amplicon_list_and_summary.RData")

amplicon <- data %>%
  tidyr::separate(amplicon_intervals, into = c("chromosome", "start", "end"), sep = "[:-]") %>%
  dplyr::filter(chromosome %in% c(1:22, "X", "Y", "M", "MT")) %>%
  dplyr::mutate(chromosome = paste0("chr", chromosome), start = as.integer(start), end = as.integer(end), idx = dplyr::row_number()) %>%
  dplyr::select(c(4:7, 1:3)) %>%
  data.table::as.data.table()

amplicon

saveRDS(amplicon[, -c("idx")], file = "data/amplicon_hg19.rds")
