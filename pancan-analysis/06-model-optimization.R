library(gcap)
library(data.table)
library(UCSCXenaShiny)

setwd(file.path(PROJ_DIR, "pancan-analysis"))

dt_class = readRDS("data/tcga_fCNA_class.rds")

proj_files = list.files("/data3/wsx_data/tcga_snp_array_gcap_result", pattern = "prediction", full.names = TRUE, all.files = TRUE)
data = purrr::map_df(proj_files, function(x) {
  message("handing", x)
  dt = readRDS(x)
  dt[sample %in% dt_class$sample]
})
#data.table::setDT(data)
uniqLen(data$sample) # 1703
dt_class = dt_class[sample %in% unique(data$sample)]
dt_type = data.table::as.data.table(tcga_clinical)[, .(sample, type)][sample %in% dt_class$sample]

dt_type[, type := factor(type, c(
  "BLCA", "BRCA", "CESC", "COAD", "DLBC", "ESCA", "GBM", "HNSC",
  "KICH", "KIRC", "KIRP", "LGG", "LIHC", "LUAD", "LUSC", "OV",
  "PRAD", "READ", "SARC", "SKCM", "STAD", "THCA", "UCEC", "UVM"
))]

data.table::setkey(data, NULL)
data = merge(data, dt_type, by = "sample")
data = mltools::one_hot(data, cols = "type")

# Different models for prediction
data$prob_xgb11 = gcap.runPrediction(data, model = readRDS("../modeling/data/xgb_v4/XGB_NF11.rds"))
data$prob_xgb32 = gcap.runPrediction(data, model = readRDS("../modeling/data/xgb_v4/XGB_NF32.rds"))
data$prob_xgb56 = gcap.runPrediction(data, model = readRDS("../modeling/data/xgb_v4/XGB_NF56.rds"))

data$prob_xgb11_d = gcap.runPrediction(data, model = readRDS("../modeling/data/xgb_v4/XGB_NF11_downsample.rds"))
data$prob_xgb32_d = gcap.runPrediction(data, model = readRDS("../modeling/data/xgb_v4/XGB_NF32_downsample.rds"))
data$prob_xgb56_d = gcap.runPrediction(data, model = readRDS("../modeling/data/xgb_v4/XGB_NF56_downsample.rds"))

data$prob_xgb32_m = gcap.runPrediction(data, model = readRDS("../modeling/data/xgb_v4/XGB_NF32_mindiff.rds"))




data_agg = data[, .(
  prob_xgb11 = max(prob_xgb11, na.rm = TRUE),
  prob_xgb32 = max(prob_xgb32, na.rm = TRUE),
  prob_xgb56 = max(prob_xgb56, na.rm = TRUE),
  prob_xgb11_d = max(prob_xgb11_d, na.rm = TRUE),
  prob_xgb32_d = max(prob_xgb32_d, na.rm = TRUE),
  prob_xgb56_d = max(prob_xgb56_d, na.rm = TRUE),
  prob_xgb32_m = max(prob_xgb32_m, na.rm = TRUE),
  CN8 = CN8[1],
  CN7_plus_CN8 = CN7[1] + CN8[1]
), by = "sample"]

dd = merge(data_agg, dt_class, by = "sample")
dd[, true_label := ifelse(class=="circular", 1L, 0L)]
dd
saveRDS(dd, file = "data/model_pred_TCGA_result.rds")

rv = data.table()
mt = c("xgb11", "xgb32", "xgb56", "xgb11_d", "xgb32_d", "xgb56_d", "xgb32_m", "CN8", "CN7_plus_CN8")
for (i in mt) {
  m = if (startsWith(i, "xgb")) paste0("prob_", i) else i
  rv = rbind(rv, data.table(
    model = i,
    auprc = get_auc(dd[[m]], dd$true_label, type = "pr")$auc.integral,
    auroc = get_auc(dd[[m]], dd$true_label, type = "roc")$auc,
    precision = ModelMetrics::precision(dd$true_label, if (max(dd[[m]]) > 1) {
      ifelse(dd[[m]] >= 1,  1, 0)
    } else dd[[m]]),
    recall = ModelMetrics::recall(dd$true_label, if (max(dd[[m]]) > 1) {
      ifelse(dd[[m]] >= 1,  1, 0)
    } else dd[[m]])
  ))
}
rv
saveRDS(rv, file = "data/sample_perf_of_model_on_TCGA_SNP_data.rds")

# XGB32 chosen
# downsample的模型差别不算大

# Determining min_prob ----------------------------------------------------
# With prob_xgb32


getf = function(i) {
  data_agg = data[, .(
    prob = max(prob_xgb32, na.rm = TRUE)
  ), by = "sample"]
  
  dd = merge(data_agg, dt_class, by = "sample")
  dd[, true_label := ifelse(class=="circular", 1L, 0L)]
  
  v = data.frame(
    acc = mean(dd$true_label == ifelse(dd$prob > i, 1, 0)),
    recall = ModelMetrics::recall(dd$true_label, dd$prob, cutoff = i),
    precision = ModelMetrics::precision(dd$true_label, dd$prob, cutoff = i),
    specificity = ModelMetrics::specificity(dd$true_label, dd$prob, cutoff = i)
  )
  print(v)
  v
}

df = purrr::map_df(
  seq(0.2, 0.95, 0.05),
  getf
)

# 关注acc和recall
df$min_prob = seq(0.2, 0.95, 0.05)
df

#          acc    recall precision specificity min_prob
# 1  0.6729301 0.9065421 0.3557457   0.6186686     0.20
# 2  0.6817381 0.9034268 0.3620474   0.6302460     0.25
# 3  0.6934821 0.9003115 0.3709884   0.6454414     0.30
# 4  0.7034645 0.8940810 0.3786280   0.6591896     0.35
# 5  0.7140341 0.8816199 0.3866120   0.6751085     0.40
# 6  0.7287140 0.8691589 0.3991416   0.6960926     0.45
# 7  0.7345860 0.8411215 0.4023845   0.7098408     0.50
# 8  0.7433940 0.8193146 0.4096573   0.7257598     0.55
# 9  0.7639460 0.7881620 0.4310051   0.7583213     0.60
# 10 0.7792132 0.7632399 0.4495413   0.7829233     0.65
# 11 0.7968291 0.7196262 0.4743326   0.8147612     0.70
# 12 0.8085731 0.6635514 0.4941995   0.8422576     0.75
# 13 0.8220787 0.5919003 0.5248619   0.8755427     0.80
# 14 0.8285379 0.4517134 0.5555556   0.9160637     0.85
# 15 0.8326483 0.3115265 0.6097561   0.9536903     0.90
# 16 0.8302995 0.1526480 0.7424242   0.9876990     0.95

saveRDS(df, file = "data/model_min_prob_perfs_xgb32.rds")

# 基本线性平缓
ggplot(df, aes(min_prob, acc)) + geom_line()
ggplot(df, aes(min_prob, specificity)) + geom_line()

# 有拐点
ggplot(df, aes(min_prob, recall)) + geom_line()
ggplot(df, aes(min_prob, precision)) + geom_line()

ggplot(df, aes(min_prob, recall)) + geom_line() +
  geom_line(aes(y = precision), color = "red")

# 基于拐点原则选择prob为0.6，这之后recall下降的很快，precision是0.9-0.95后才有明显增加
# 也就是默认0.6，如果强调ecDNA的准确性使用0.95作为阈值

# 看下xgb11，结果差不多
getf = function(i) {
  data_agg = data[, .(
    prob = max(prob_xgb11, na.rm = TRUE)
  ), by = "sample"]
  
  dd = merge(data_agg, dt_class, by = "sample")
  dd[, true_label := ifelse(class=="circular", 1L, 0L)]
  
  v = data.frame(
    acc = mean(dd$true_label == ifelse(dd$prob > i, 1, 0)),
    recall = ModelMetrics::recall(dd$true_label, dd$prob, cutoff = i),
    precision = ModelMetrics::precision(dd$true_label, dd$prob, cutoff = i),
    specificity = ModelMetrics::specificity(dd$true_label, dd$prob, cutoff = i)
  )
  print(v)
  v
}

df = purrr::map_df(
  seq(0.2, 0.95, 0.05),
  getf
)

# 关注acc和recall
df$min_prob = seq(0.2, 0.95, 0.05)
df


# acc    recall precision specificity min_prob
# 1  0.6617733 0.9127726 0.3483948   0.6034732     0.20
# 2  0.6688197 0.9127726 0.3534379   0.6121563     0.25
# 3  0.6776277 0.9127726 0.3599509   0.6230101     0.30
# 4  0.6846741 0.9003115 0.3639798   0.6345876     0.35
# 5  0.6928949 0.8971963 0.3701799   0.6454414     0.40
# 6  0.6993541 0.8847352 0.3741765   0.6562952     0.45
# 7  0.7052261 0.8629283 0.3768707   0.6685962     0.50
# 8  0.7234292 0.8566978 0.3928571   0.6924747     0.55
# 9  0.7428068 0.8348910 0.4104135   0.7214182     0.60
# 10 0.7563124 0.7912773 0.4219269   0.7481910     0.65
# 11 0.7692308 0.7196262 0.4325843   0.7807525     0.70
# 12 0.7903699 0.6386293 0.4596413   0.8256151     0.75
# 13 0.8073987 0.5420561 0.4901408   0.8690304     0.80
# 14 0.8209043 0.4330218 0.5305344   0.9109986     0.85
# 15 0.8256019 0.2897196 0.5740741   0.9500724     0.90
# 16 0.8291251 0.1401869 0.7500000   0.9891462     0.95

saveRDS(df, file = "data/model_min_prob_perfs_xgb11.rds")

# 基本线性平缓
ggplot(df, aes(min_prob, acc)) + geom_line()
ggplot(df, aes(min_prob, specificity)) + geom_line()

# 有拐点
ggplot(df, aes(min_prob, recall)) + geom_line()
ggplot(df, aes(min_prob, precision)) + geom_line()

ggplot(df, aes(min_prob, recall)) + geom_line() +
  geom_line(aes(y = precision), color = "red")

# 选择prob为0.6，这之后recall下降的很快，precision是0.9-0.95后才有明显增加
# 也就是默认0.6，如果强调ecDNA的准确性使用0.95作为阈值


