library(gcap)
library(data.table)

setwd(file.path(PROJ_DIR, "pancan-analysis"))

dt = readRDS("data/pcawg_copynumber_sp.rds")
colnames(dt) = c("sample", "chromosome", "start", "end", "total_cn", "minor_cn")
dt[, chromosome := paste0("chr", chromosome)]

cli = readRDS("data/pcawg_sample_info.rds")
cli2 = cli[, .(sample, purity, ploidy, age = donor_age_at_diagnosis, gender = donor_sex, type = cancer_type)]
cli2[, gender := ifelse(gender == "female", "XX", "XY")]

dt2 = merge(dt, cli2, by = "sample", all.x = TRUE)

rm(dt, cli, cli2); gc()
saveRDS(dt2, file = "data/pcawg_ascat.rds")

library(gcap)
library(data.table)
setwd(file.path(PROJ_DIR, "pancan-analysis"))

dt2 = readRDS("data/pcawg_ascat.rds")

#dt3 = copy(dt2)
#dt2 = dt2[type == "Biliary-AdenoCA"]
dt_list = lapply(split(dt2, dt2$type), function(data) {
  #data.table::setDTthreads(threads = 1) # https://github.com/HenrikBengtsson/future/issues/343
  type = data$type[1]
  #model = "XGB11"
  model = readRDS("../modeling/data/xgb_v3/XGB_NF11.rds")
  data$type = NULL
  message("Handling type ", type, " with sample size ", length(unique(data$sample)))
  
  rv = suppressMessages(
    gcap.ASCNworkflow(data, 
                      outdir = "/data3/wsx_data/pcawg_gcap_result2",
                      result_file_prefix = paste0("PCAWG_", make.names(type)),
                      genome_build = "hg19",
                      model = model, 
                      tightness = 1, gap_cn = 3)
  )
})
