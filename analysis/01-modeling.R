PROJ_DIR = "~/gcap-analysis/manuscript/"
setwd(PROJ_DIR)


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# CNV source platform comparison
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

comp_df_all <- data.table::rbindlist(
  list(
    wgs_snp = readRDS("../preprocessing-and-EDA/data/cnv_comparison_result_between_TCGA_WGS_and_SNParray_allgenes.rds"),
    wgs_wes = readRDS("../preprocessing-and-EDA/data/cnv_comparison_result_between_TCGA_WGS_and_WES_allgenes.rds"),
    snp_wes = readRDS("../preprocessing-and-EDA/data/cnv_comparison_result_between_TCGA_WES_and_SNParray_allgenes.rds")
  ),
  idcol = "platform_pair"
)

library(ggblanket)

plot_comparison <- function(comp_df) {
  count_df <- comp_df %>% dplyr::count(platform_pair)
  pair_levels <- c("wgs_snp", "wgs_wes", "snp_wes")
  names(pair_levels) <- c(
    paste0("WGS and SNParray", paste0("\n(n=", count_df %>% dplyr::filter(platform_pair == "wgs_snp") %>% dplyr::pull(n), ")")),
    paste0("WGS and WES", paste0("\n(n=", count_df %>% dplyr::filter(platform_pair == "wgs_wes") %>% dplyr::pull(n), ")")),
    paste0("WES and SNParray", paste0("\n(n=", count_df %>% dplyr::filter(platform_pair == "snp_wes") %>% dplyr::pull(n), ")"))
  )
  
  df <- comp_df %>%
    tidyr::pivot_longer(cols = c("sim_total", "rmse_total", "sim_minor", "rmse_minor"), names_to = "measure", values_to = "score") %>%
    dplyr::mutate(
      platform_pair = forcats::fct_recode(platform_pair, !!!pair_levels),
      measure = forcats::fct_recode(measure,
                                    `Total CN similarity` = "sim_total",
                                    `Total CN RMSE` = "rmse_total",
                                    `Minor CN similarity` = "sim_minor",
                                    `Minor CN RMSE` = "rmse_minor"
      )
    )
  
  library(ggplot2)
  p1 = gg_boxplot(df %>% dplyr::filter(measure %in% c("Total CN RMSE", "Minor CN RMSE")), 
                  platform_pair, score, col = platform_pair, y_expand = c(0, 0)) +
    facet_wrap(~measure, scales = "free", nrow = 1) +
    labs(x = NULL, y = NULL)
  p2 = gg_boxplot(df %>% dplyr::filter(measure %in% c("Total CN similarity", "Minor CN similarity")), 
                  platform_pair, score, col = platform_pair, y_expand = c(0, 0)) +
    facet_wrap(~measure, scales = "free", nrow = 1) +
    labs(x = NULL, y = NULL)
  cowplot::plot_grid(p1, p2, nrow = 2)
}

pdf("plots/cnv_platform_comparison.pdf", width = 10, height = 5)
plot_comparison(comp_df_all)
dev.off()

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Modeling metrics and comparison
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# TOP 20的表格可以输出

top20 = list(
  xgb11 = readRDS("../modeling/data/xgb_v4/XGB_NF11_perf_top20.rds"),
  xgb32 = readRDS("../modeling/data/xgb_v4/XGB_NF32_perf_top20.rds"),
  xgb56 = readRDS("../modeling/data/xgb_v4/XGB_NF56_perf_top20.rds"),
  # downsample
  
  xgb11_downsample = readRDS("../modeling/data/xgb_v4/XGB_NF11_downsample_perf_top20.rds"),
  xgb32_downsample = readRDS("../modeling/data/xgb_v4/XGB_NF32_downsample_perf_top20.rds"),
  xgb56_downsample = readRDS("../modeling/data/xgb_v4/XGB_NF56_downsample_perf_top20.rds")
)

top_aucpr = purrr::map2_df(top20, names(top20),
                           ~data.frame(model = .y, 
                                       auPRC_mean = .$test_aucpr_mean[1], 
                                       auPRC_sd = .$test_aucpr_sd[1]))
top_aucpr
gg_bar(top_aucpr, x = model, y = auPRC_mean, col = model, stat = "identity")

p = gg_blank(
  top_aucpr %>% 
    dplyr::mutate(
      auPRC_mean = round(auPRC_mean, 3),
      auPRC_sd = round(auPRC_sd, 3),
      type = ifelse(endsWith(model, "downsample"), "downsample", "full"),
      model = sub("_downsample", "", model)
    ), 
  x = model, y = auPRC_mean, col = type, pal = RColorBrewer::brewer.pal(3, "Dark2"),
  facet = type, facet_scales = "free",
  label = auPRC_mean,
  ymin = auPRC_mean,
  ymax = auPRC_mean + auPRC_sd,
  yend = auPRC_mean + auPRC_sd + 0.05,
  y_include = 0,
  y_title = "auPRC (area under precision-recall curve)"
) +
  geom_col(width = 0.75, alpha = 0.9) +
  geom_errorbar(width = 0.1, colour = pal_na()) +
  geom_text(aes(y = auPRC_mean + auPRC_sd + 0.02), fill = NA, size = 3)
p
ggsave("plots/model_cv_performance.pdf", p, width = 8, height = 4)

# downsample模型使用的数据并不一样
# 统一使用全部的数据
get_perf <- function(bst, idx) {
  old_wd = getwd()
  wd = "/data3/wsx_data/gcap-analysis"
  setwd(wd)
  on.exit(setwd(old_wd))
  ## Use the idx to determine the target
  if (idx == 1L) {
    rm_type = TRUE
    rm_cli = FALSE
    rm_cns = FALSE
  } else if (idx == 2L) {
    rm_type = TRUE
    rm_cli = TRUE
    rm_cns = TRUE
  } else if (idx == 3L) {
    rm_type = FALSE
    rm_cli = FALSE
    rm_cns = FALSE
  }
  
  library(stats)
  suppressMessages(library(xgboost))
  stopifnot(requireNamespace("gcap", quietly = TRUE))
  
  ## Preprocess data
  source("lib/get_data.R")
  data = get_model_data2(rm_type = rm_type, rm_cli = rm_cli, rm_cns = rm_cns,
                              version = "v4", fully_random = TRUE)
  dtrain = xgb.DMatrix(data$data[, !colnames(data$data) %in% "y"],
                       label = data$data[, "y"])
  
  y_pred = predict(bst, dtrain)
  y_true = xgboost::getinfo(dtrain, "label")
  rv = data.frame(
    auPRC = gcap::get_auc(y_pred, y_true)$auc.integral,
    auROC = gcap::get_auc(y_pred, y_true, type = "roc")$auc,
    precision = ModelMetrics::precision(y_true, y_pred),
    recall = ModelMetrics::recall(y_true, y_pred)
  )
  print(rv)
  rv
}


ml = list(
  xgb11 = readRDS("../modeling/data/xgb_v4/XGB_NF11.rds"),
  xgb32 = readRDS("../modeling/data/xgb_v4/XGB_NF32.rds"),
  xgb56 = readRDS("../modeling/data/xgb_v4/XGB_NF56.rds"),
  # downsample
  xgb11_downsample = readRDS("../modeling/data/xgb_v4/XGB_NF11_downsample.rds"),
  xgb32_downsample = readRDS("../modeling/data/xgb_v4/XGB_NF32_downsample.rds"),
  xgb56_downsample = readRDS("../modeling/data/xgb_v4/XGB_NF56_downsample.rds")
)

perf = purrr::map2_df(ml, c(2, 1, 3, 2, 1, 3), get_perf)
perf$model = names(ml)
perf

dir.create("data")
saveRDS(perf, file = "data/model_perf_on_full_modeling_data.rds")

# Pick XGB11 model
eval_log = readRDS(file.path("../modeling/", top20$xgb11$file[1]))[[1]]
eval_log$evaluation_log

eval_df = rbind(
  data.table(type = "train", iter = eval_log$evaluation_log$iter,
             aucpr_mean = eval_log$evaluation_log$train_aucpr_mean,
             aucpr_sd = eval_log$evaluation_log$train_aucpr_std),
  data.table(type = "test", iter = eval_log$evaluation_log$iter,
             aucpr_mean = eval_log$evaluation_log$test_aucpr_mean,
             aucpr_sd = eval_log$evaluation_log$test_aucpr_std)
)
eval_df[, type := factor(type, c("train", "test"))]

p = gg_blank(
  eval_df,
  x = iter, y = aucpr_mean, col = type, pal = RColorBrewer::brewer.pal(3, "Dark2"),
  ymin = aucpr_mean - aucpr_sd,
  ymax = aucpr_mean + aucpr_sd,
  yend = aucpr_mean + aucpr_sd,
  y_include = 0,
  y_limits = c(0, 1),
  col_legend_place = "t",
  y_title = "auPRC (area under precision-recall curve)",
  x_title = "Iteration"
) +
  geom_line(width = 0.75, alpha = 0.9) +
  geom_errorbar(width = 0.1, colour = pal_na()) +
  geom_vline(xintercept = eval_log$best_ntreelimit, linetype = 2)
p + cowplot::theme_cowplot()

ggsave("plots/XGB11_cv_rounds2.pdf", p, width = 6, height = 4)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Sample level classifier
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
library(gcap)
library(data.table)
source("../lib/get_data.R")
all_samples = get_model_data2(only_samples = TRUE)

data_label = readxl::read_excel("../preprocessing-and-EDA/data/ICGC-ecDNA.xlsx", sheet = 3)
data_label = data_label %>% dplyr::filter(sample_barcode %in% all_samples) %>% 
  dplyr::select(sample = sample_barcode, y_label = sample_classification) %>% 
  data.table::as.data.table()

data_label[, y_true := fcase(y_label %in% c("Circular"), 1L, default = 0L)]

data = readRDS("../modeling/data/train_data_v4_raw.rds")
colnames(data)[3] <- "total_cn"
data$prob = gcap.runPrediction(data, model = "XGB11")
# data2 = gcap.runScoring(data)
# data2 = data2$fCNA$sample_summary
# df = merge(data_label, data2[, .(sample, y_pred = ifelse(class == "circular", 1, 0))], by = "sample")

# Try 0.5, 0.6, 0.7, 0.8, 0.9 and 0.6 got the best score
df = merge(data_label, data[, .(y_pred = ifelse(max(prob, na.rm = TRUE) > 0.6, 1, 0)), by = "sample"], by = "sample")
df
get_auc(df$y_pred, df$y_true, type = "roc")$auc


sper = data.table(
  measure = c("auPRC", "auROC", "sensitivity", "precision", "specificity"),
  score = c(
    get_auc(df$y_pred, df$y_true, type = "pr")$auc.integral,
    ModelMetrics::auc(df$y_true, df$y_pred),
    ModelMetrics::sensitivity(df$y_true, df$y_pred),
    ModelMetrics::precision(df$y_true, df$y_pred),
    ModelMetrics::specificity(df$y_true, df$y_pred)
  )
)
sper

saveRDS(sper, file = "data/sample_class_perf_on_modeling_data.rds")

p = sper %>% 
  ggplot(aes(x = measure, y = score)) + 
  geom_col(fill = "steelblue") +
  geom_text(aes(label = round(score, 3)), nudge_y = 0.05) +
  cowplot::theme_cowplot() +
  xlab(NULL) + ylab("Performance score")
p

ggsave("plots/model_sample_performance.pdf", p, width = 6, height = 4)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Models with different folds
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

rv_df = list(
  XGB11 = readRDS("data/NF11_with_diff_folds.rds"),
  XGB32 = readRDS("data/NF32_with_diff_folds.rds"),
  XGB56 = readRDS("data/NF56_with_diff_folds.rds")
)

rv_lg = readRDS("data/model_perf_logistic.rds")

rv_list = lapply(rv_df, function(z) {
  dd = purrr::map_df(z, function(x) {
    rv = purrr::map_df(x, ~.$evaluation_log[nrow(.$evaluation_log), ], .id = "repeat")
    print(rv)
    rv
  }, .id = "folds")
  
  dd[, aucpr := round(test_aucpr_mean, 3)]
  dd
})

rv_list$XGB11

library(ggpubr)

plist = purrr::map(rv_list,
           ~ggpubr::ggboxplot(.,
                              x = "folds", y = "aucpr",
                              fill = "steelblue",
                              width = 0.5) +
             labs(y = "auPRC"))



ggsave("plots/auPRC_with_different_folds_XGB11.pdf", plot = plist[[1]], width = 5, height = 3)
ggsave("plots/auPRC_with_different_folds_XGB32.pdf", plot = plist[[2]], width = 5, height = 3)
ggsave("plots/auPRC_with_different_folds_XGB56.pdf", plot = plist[[3]], width = 5, height = 3)


rv_lg2 = purrr::map_df(rv_lg, function(x) {
  rv = purrr::map_df(x, ~data.table::data.table(aucpr = mean(.),
                                                aucpr_sd = sd(.)), .id = "repeat")
  print(rv)
  rv
}, .id = "folds")
rv_lg2[, folds := as.integer(folds)]
rv_lg2[, model := "logistic"]

rv_data = rbindlist(rv_list, idcol = "model")
rv_data

rv_data2 = rbind(rv_data[, c("model", "folds", "repeat", "aucpr")],
                 rv_lg2[, c("model", "folds", "repeat", "aucpr")])

# Add logistic result

p = ggpubr::ggboxplot(rv_data2,
                  x = "folds", y = "aucpr",
                  fill = "model",
                  width = 0.5) +
  labs(y = "auPRC")
p

ggsave("plots/auPRC_with_different_folds_XGB_all.pdf", plot = p, width = 7, height = 3)
