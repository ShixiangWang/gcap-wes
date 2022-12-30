library(IDConverter)
library(dplyr)

setwd(PROJ_DIR)

data_label = readxl::read_excel("preprocessing-and-EDA/data/ICGC-ecDNA.xlsx", sheet = 3)
data_label = data_label %>% 
  #dplyr::filter(startsWith(sample_barcode, "SA")) %>%
  dplyr::select(sample = sample_barcode, label = sample_classification, lineage) %>% 
  data.table::as.data.table()

data_label

options(IDConverter.datapath = system.file("extdata", package = "IDConverter"))

data_label$sample_sp = convert_icgc(data_label$sample, "icgc_sample_id", "icgc_specimen_id")
sum(is.na(data_label$sample_sp))

data_label

library(data.table)
icgc2 = load_data("icgc")[, list(submitted_sample_id, icgc_specimen_id)]
icgc2[, submitted_sample_id := substr(submitted_sample_id, 1, 15)]
icgc2 = unique(icgc2)

data_label[, sample_sp := ifelse(is.na(sample_sp),
                                convert_custom(sample, "submitted_sample_id", "icgc_specimen_id", dt = icgc2),
                                sample_sp)]
data_label = data_label[!is.na(sample_sp)]

data_label[, class := fcase(label %in% c("Circular"), "circular",
                            label %in% c("No-fSCNA"), "nofocal",
                            default = "noncircular")]

saveRDS(data_label, file = "pancan-analysis/data/pcawg_fCNA_class.rds")

# 看下肠癌的circular情况
dput(data_label[lineage == "Colorectal" & class == "circular"]$sample_sp)

icgc2[icgc_specimen_id %in% c("SP17953", "SP20419", "SP18430", "SP17581")]
# submitted_sample_id icgc_specimen_id
# 1:     TCGA-AZ-6601-01          SP17581
# 2:     TCGA-A6-3810-01          SP20419
# 3:     TCGA-A6-5656-01          SP18430
# 4:     TCGA-A6-2677-01          SP17953

# PCAWG的拷贝数数据没有包含它们

