setwd(file.path(PROJ_DIR, "preprocessing-and-EDA/"))

library(readxl)
library(tidyverse)

data <- read_excel("data/ICGC-ecDNA.xlsx", sheet = 1)
all_tcga_tumor_info <- read_excel("data/ICGC-ecDNA.xlsx", sheet = 3) %>%
  dplyr::filter(startsWith(sample_barcode, "TCGA"), tumor_or_normal == "tumor") %>%
  dplyr::select(sample_barcode, sample_classification) %>%
  unique()

data <- separate_rows(data, amplicon_intervals, sep = ",")

# Total used sample number is 3212, as described in Nature Genetics (2020) paper
length(table(data$sample_barcode)) # 1348/3212

data_summary <- data %>%
  group_by(sample_barcode, amplicon_classification) %>%
  summarise(N = n()) %>%
  pivot_wider(sample_barcode, names_from = "amplicon_classification", values_from = "N", values_fill = 0) %>%
  arrange(desc(Circular))

data_summary_tcga <- data_summary %>%
  filter(startsWith(sample_barcode, "TCGA")) %>%
  arrange(desc(Circular))
# This need appending samples with amplicons available
left_ids <- unique(setdiff(all_tcga_tumor_info$sample_barcode, data_summary_tcga$sample_barcode))
data_summary_tcga <- dplyr::bind_rows(data_summary_tcga, dplyr::tibble(id = left_ids, a = 0L, b = 0L, c = 0L, d = 0L) %>%
  setNames(colnames(data_summary_tcga)))

data_summary_with_amplicon <- data_summary

save(data, data_summary_with_amplicon, data_summary_tcga, all_tcga_tumor_info, file = "data/pancan_amplicon_list_and_summary.RData")

# Select samples for downstream analysis ---------
load("data/pancan_amplicon_list_and_summary.RData")
length(data_summary_tcga$sample_barcode[data_summary_tcga$Circular > 3])
# Select samples with at least 4 circular amplicon regions
paste(substr(data_summary_tcga$sample_barcode[data_summary_tcga$Circular > 3], 1, 12), collapse = ",")

# Add non-ecDNA samples
data2 = data_summary_tcga %>% filter(Circular == 0) %>% 
  left_join(all_tcga_tumor_info) %>% 
  ungroup()
paste(substr(data2$sample_barcode, 1, 12), collapse = ",")

info = readr::read_tsv("../gdc-manifest/gdc_sample_sheet.2021-12-23.tsv")

data2 = data2 %>% filter(sample_barcode %in% substr(info$`Sample ID`, 1, 15))
data2 = data2 %>% filter(sample_classification %in% c("No-fSCNA", "Linear"))


table(data2$sample_classification )

set.seed(2021)
data3 = sample_n(data2 %>% filter(sample_classification == "No-fSCNA"), 30)
data4 = sample_n(data2 %>% filter(sample_classification == "Linear"), 30)

data5 = bind_rows(data3, data4)
paste(substr(data5$sample_barcode, 1, 12), collapse = ",")


# Select samples with 1-3 circular amplicon counts

paste(substr(data_summary_tcga$sample_barcode[data_summary_tcga$Circular > 0 & data_summary_tcga$Circular <= 3], 1, 12), collapse = ",")

# 2021-10-15
