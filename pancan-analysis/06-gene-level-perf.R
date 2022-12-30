library(gcap)
library(data.table)

setwd(file.path(PROJ_DIR, "pancan-analysis"))

dt_class = readRDS("data/tcga_fCNA_class.rds")

proj_files = list.files("/data3/wsx_data/tcga_snp_array_gcap_result", pattern = "prediction", full.names = TRUE, all.files = TRUE)
data = purrr::map_df(proj_files, function(x) {
  message("handing", x)
  dt = readRDS(x)
  dt[sample %in% dt_class$sample]
})
#data.table::setDT(data)
uniqLen(data$sample)
dt_class = dt_class[sample %in% unique(data$sample)]

library(UCSCXenaShiny)
dt_type = data.table::as.data.table(tcga_clinical)[, .(sample, type)][sample %in% dt_class$sample]

dt_type[, type := factor(type, c(
  "BLCA", "BRCA", "CESC", "COAD", "DLBC", "ESCA", "GBM", "HNSC",
  "KICH", "KIRC", "KIRP", "LGG", "LIHC", "LUAD", "LUSC", "OV",
  "PRAD", "READ", "SARC", "SKCM", "STAD", "THCA", "UCEC", "UVM"
))]

data.table::setkey(data, NULL)
data = merge(data, dt_type, by = "sample")
data = mltools::one_hot(data, cols = "type")

# XGB
data$prob_xgb11 = gcap.runPrediction(data, model = readRDS("../modeling/data/xgb_v4/XGB_NF11.rds"))
data$prob_xgb32 = gcap.runPrediction(data, model = readRDS("../modeling/data/xgb_v4/XGB_NF32.rds"))
data$prob_xgb56 = gcap.runPrediction(data, model = readRDS("../modeling/data/xgb_v4/XGB_NF56.rds"))

data$prob_xgb11_d = gcap.runPrediction(data, model = readRDS("../modeling/data/xgb_v4/XGB_NF11_downsample.rds"))
data$prob_xgb32_d = gcap.runPrediction(data, model = readRDS("../modeling/data/xgb_v4/XGB_NF32_downsample.rds"))
data$prob_xgb56_d = gcap.runPrediction(data, model = readRDS("../modeling/data/xgb_v4/XGB_NF56_downsample.rds"))

ascat = readRDS("data/tcga_ascat_tidy.rds")
ascat = ascat[, .(chr = paste0("chr", chromosome), start, end, total_cn, minor_cn, sample, CNsignatureMapping)]
setkey(ascat, NULL)
regions = modules::use("../lib/regions.R")
# ref = readRDS(
#   system.file(
#     "extdata", paste0("hg38", "_target_genes.rds"),
#     package = "gcap", mustWork = TRUE
#   )
# )
ref_hg19 = readRDS(
  system.file(
    "extdata", paste0("hg19", "_target_genes.rds"),
    package = "gcap", mustWork = TRUE
  )
)
ascat_sigmap = regions$overlaps(ascat, ref_hg19)[intersect_ratio >= 1]
ascat_sigmap

dd = merge(data[, .(sample, gene_id, total_cn, minor_cn, ploidy, 
                    prob_xgb11, prob_xgb32, prob_xgb56,
                    prob_xgb11_d, prob_xgb32_d, prob_xgb56_d)],
           ascat_sigmap[, .(sample, gene_id, CNsignatureMapping)], 
           by = c("sample", "gene_id"))

amplicon <- readRDS("../preprocessing-and-EDA/data/amplicon_hg19.rds")
amplicon = amplicon[sample_barcode %in% dt_class$sample]
uniqLen(amplicon$sample_barcode)
amplicon_dt = regions$overlaps(amplicon, ref_hg19)[intersect_ratio >= 1]

dd = merge(dd[gene_id %in% ref_hg19$gene_id],
           amplicon_dt[, .(sample = sample_barcode, gene_id, amplicon_classification)], 
           by = c("sample", "gene_id"), all.x = TRUE)

# Remove unreasonable records
dd2 = dd[!(((total_cn < ploidy + 4) & !is.na(amplicon_classification)) | (is.na(amplicon_classification) & (total_cn >= ploidy + 4)))]

uniqLen(dd2$sample)
# [1] 1703
nrow(dd2)
# [1] 32354524

dd2[, label := ifelse(amplicon_classification %in% "Circular", 1L, 0L)]
table(dd2$label)
dd2[, CN8 := ifelse(CNsignatureMapping %in% "CN8", 1L, 0L)]
dd2[, CN7_plus_CN8 := ifelse(CNsignatureMapping %in% c("CN7", "CN8"), 1L, 0L)]


rv = data.table()
mt = c("xgb11", "xgb32", "xgb56", "xgb11_d", "xgb32_d", "xgb56_d", "CN8", "CN7_plus_CN8")
for (i in mt) {
  m = if (startsWith(i, "xgb")) paste0("prob_", i) else i
  rv = rbind(rv, data.table(
    model = i,
    auprc = get_auc(dd2[[m]], dd2$label, type = "pr")$auc.integral,
    auroc = get_auc(dd2[[m]], dd2$label, type = "roc")$auc,
    precision = ModelMetrics::precision(dd2$label, dd2[[m]]),
    recall = ModelMetrics::recall(dd2$label, dd2[[m]])
  ))
}
rv

saveRDS(rv, file = "data/gene_perf_of_model_on_TCGA_SNP_data.rds")
