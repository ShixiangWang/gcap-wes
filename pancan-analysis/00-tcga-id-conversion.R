library(IDConverter)
library(dplyr)

setwd(PROJ_DIR)

data_label = readxl::read_excel("preprocessing-and-EDA/data/ICGC-ecDNA.xlsx", sheet = 3)
data_label = data_label %>% 
  dplyr::filter(startsWith(sample_barcode, "TCGA")) %>%
  dplyr::filter(as.integer(substr(sample_barcode, 14, 15)) < 10) %>% 
  dplyr::select(sample = sample_barcode, label = sample_classification) %>% 
  data.table::as.data.table()

data_label

library(data.table)
data_label[, class := fcase(label %in% c("Circular"), "circular",
                            label %in% c("No-fSCNA"), "nofocal",
                            default = "noncircular")]

saveRDS(data_label, file = "pancan-analysis/data/tcga_fCNA_class.rds")


