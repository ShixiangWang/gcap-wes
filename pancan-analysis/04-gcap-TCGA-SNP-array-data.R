library(gcap)
library(data.table)

setwd(file.path(PROJ_DIR, "pancan-analysis"))

dt <- readRDS("../preprocessing-and-EDA/data/tcga_ascat.rds")
data.table::setnames(dt, "cancer_type", "type")
data.table::setnames(dt, "chr", "chromosome")
data.table::setnames(dt, "sex", "gender")

library(UCSCXenaShiny)
library(dplyr)
cli = tcga_clinical %>% 
  select(sample, age = age_at_initial_pathologic_diagnosis) %>%
  data.table::as.data.table() %>% 
  unique()

dt2 = merge(dt, cli, by = "sample", all.x = TRUE)
rm(dt); gc()

saveRDS(dt2, file = "data/tcga_ascat_tidy.rds")

dt2 = readRDS("data/tcga_ascat_tidy.rds")

dt_list = lapply(split(dt2, dt2$type), function(data) {
  type = data$type[1]
  #model = "XGB11"
  model = readRDS("../modeling/data/xgb_v3/XGB_NF11.rds")
  data$type = NULL
  message("Handling type ", type, " with sample size ", length(unique(data$sample)))
  Sys.sleep(1)

  rv = gcap.ASCNworkflow(data, 
                    outdir = "/data3/wsx_data/tcga_snp_array_gcap_result2",
                    result_file_prefix = paste0("TCGA_SNP_", type),
                    genome_build = "hg19",
                    model = model,
                    tightness = 1,
                    gap_cn = 3)
})

