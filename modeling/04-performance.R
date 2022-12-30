PROJ_DIR = "/data3/wsx_data/gcap-analysis/modeling/"
setwd(PROJ_DIR)

suppressMessages(source("../lib/performance.R"))

model_list = list(
  NF32 = readRDS("data/xgb_v3/XGB_NF32.rds"),
  NF11 = readRDS("data/xgb_v3/XGB_NF11.rds"),
  NF56 = readRDS("data/xgb_v3/XGB_NF56.rds"),
  NF32_subset = readRDS("data/xgb_v3/XGB_NF32_subset.rds"),
  NF32_downsample = readRDS("data/xgb_v3/XGB_NF32_downsample.rds")
)

perf = list()

perf$NF32 = get_perf(model_list$NF32, 1)
perf$NF11 = get_perf(model_list$NF11, 2)
perf$NF56 = get_perf(model_list$NF56, 3)
perf$NF32_subset = get_perf(model_list$NF32_subset, 1)
perf$NF32_downsample = get_perf(model_list$NF32_downsample, 1)

# Save the result for further use
saveRDS(perf, file = "data/xgb_v3/performance_of_model_list.rds")

# Check the data
str(perf, max.level = 2)

str(perf$NF32$perf_overall, max.level = 3)

library(ggplot2)

perf = readRDS("data/xgb_v3/performance_of_model_list.rds")

p = plot_perf(perf$NF32$perf_overall)
p
ggsave("plots/xgb32_ROC_curve.pdf", p, width = 5, height = 4)

p = plot_perf(perf$NF32$perf_overall, "auPR")
p
ggsave("plots/xgb32_PR_curve.pdf", p, width = 5, height = 4)

p1 = plot_perf(perf$NF11$perf_overall)
p2 = plot_perf(perf$NF11$perf_overall, "auPR")
p1
p2

ggsave("plots/xgb11_ROC_curve.pdf", p1, width = 5, height = 4)
ggsave("plots/xgb11_PR_curve.pdf", p2, width = 5, height = 4)

p1 = plot_perf(perf$NF56$perf_overall)
p2 = plot_perf(perf$NF56$perf_overall, "auPR")

p1
p2

ggsave("plots/xgb56_ROC_curve.pdf", p1, width = 5, height = 4)
ggsave("plots/xgb56_PR_curve.pdf", p2, width = 5, height = 4)


p1 = plot_perf(perf$NF32_subset$perf_overall)
p2 = plot_perf(perf$NF32_subset$perf_overall, "auPR")

p1
p2

ggsave("plots/xgb32_subset_ROC_curve.pdf", p1, width = 5, height = 4)
ggsave("plots/xgb32_subset_PR_curve.pdf", p2, width = 5, height = 4)

p1 = plot_perf(perf$NF32_downsample$perf_overall)
p2 = plot_perf(perf$NF32_downsample$perf_overall, "auPR")

p1
p2

ggsave("plots/xgb32_downsample_ROC_curve.pdf", p1, width = 5, height = 4)
ggsave("plots/xgb32_downsample_PR_curve.pdf", p2, width = 5, height = 4)


NF11 = c(perf$NF11$perf_overall$train$PR$auc, perf$NF11$perf_overall$test$PR$auc)
NF32 = c(perf$NF32$perf_overall$train$PR$auc, perf$NF32$perf_overall$test$PR$auc)
NF56 = c(perf$NF56$perf_overall$train$PR$auc, perf$NF56$perf_overall$test$PR$auc)

data = dplyr::tibble(
  set = rep(c("train", "test"), 3) %>% forcats::fct_inorder(),
  nfeature = c(rep(11, 2), rep(32, 2), rep(56, 2)),
  aucpr = c(NF11, NF32, NF56)
)


p = ggplot(data, aes(y = aucpr, x = factor(nfeature))) +
  geom_point(aes(shape = set, color = set), size = 3) +
  cowplot::theme_cowplot() +
  labs(x= "nfeature", y = "AUCPR")
p

ggsave("plots/xgb_comparison_dotplot.pdf", p, width = 5, height = 4)

library(forcats)

setplot(12, 4)
p = ggplot(perf$NF11$perf_stats %>% 
             dplyr::mutate(data_type = fct_inorder(data_type)),
           aes(x = cutoff, y = score, color = data_type, shape = data_type)) +
  geom_point(alpha = 0.6) +
  facet_wrap(~measure, nrow = 1, scales = "free") +
  cowplot::theme_cowplot() +
  scale_color_manual(values = c(train = "red", test = "blue")) +
  labs(x = "prob threshold", color = "data set", shape = "data set") +
  ylim(0, 1)
p

cowplot::save_plot("plots/performance_stats_under_different_prob_threshold_NF11.pdf", plot = p, base_asp = 1/3, base_width = 12)

setplot(12, 4)
p = ggplot(perf$NF32$perf_stats %>% 
             dplyr::mutate(data_type = fct_inorder(data_type)),
           aes(x = cutoff, y = score, color = data_type, shape = data_type)) +
  geom_point(alpha = 0.6) +
  facet_wrap(~measure, nrow = 1, scales = "free") +
  cowplot::theme_cowplot() +
  scale_color_manual(values = c(train = "red", test = "blue")) +
  labs(x = "prob threshold", color = "data set", shape = "data set") +
  ylim(0, 1)
p

cowplot::save_plot("plots/performance_stats_under_different_prob_threshold_NF32.pdf", plot = p, base_asp = 1/3, base_width = 12)

setplot(12, 4)
p = ggplot(perf$NF56$perf_stats %>% 
             dplyr::mutate(data_type = fct_inorder(data_type)),
           aes(x = cutoff, y = score, color = data_type, shape = data_type)) +
  geom_point(alpha = 0.6) +
  facet_wrap(~measure, nrow = 1, scales = "free") +
  cowplot::theme_cowplot() +
  scale_color_manual(values = c(train = "red", test = "blue")) +
  labs(x = "prob threshold", color = "data set", shape = "data set") +
  ylim(0, 1)
p

library(xgboost)

plot_xgb = function(fit) {
  xgb.ggplot.importance(xgb.importance(model = fit))
}


p1 = plot_xgb(model_list$NF11)
p2 = plot_xgb(model_list$NF32)
p3 = plot_xgb(model_list$NF56)

ggsave("plots/xgb11_importance.pdf", p1, width = 3, height = 4)
ggsave("plots/xgb32_importance.pdf", p2, width = 4, height = 5)
ggsave("plots/xgb56_importance.pdf", p3, width = 4, height = 6)

library(xgboost)
library(ggplot2)
ggplot(xgb.importance(
  model = readRDS("../modeling/data/xgb_v3/XGB_NF11.rds")) %>% 
    dplyr::mutate(Feature = factor(Feature, Feature))) +
  geom_bar(aes(x = Feature, y = Gain), stat = "identity", fill = "blue") +
  cowplot::theme_cowplot() +
  ggplot2::scale_fill_manual(values = c("blue", "blue", "blue")) +
  ggplot2::theme(legend.position = "none")  +
  labs(x = NULL, y = "Importance") +
  ggpubr::rotate_x_text(45) -> p
p
ggsave("plots/xgb11_importance2.pdf", p, width = 6, height = 4)


p = ggplot(perf$NF11$perf_stats %>% 
             dplyr::filter(measure != "F1") %>%
             dplyr::mutate(data_type = fct_inorder(data_type)),
           aes(x = cutoff, y = score, color = data_type, shape = data_type)) +
  geom_point(alpha = 0.6) +
  facet_wrap(~measure, nrow = 1, scales = "free") +
  cowplot::theme_cowplot() +
  scale_color_manual(values = c(train = "red", test = "blue")) +
  labs(x = "Probability", y = "Performance score",
       color = "data set", shape = "data set") +
  ylim(0, 1) +
  theme(legend.position = "top")
p

cowplot::save_plot("../manuscript/plots/gene_performance_score.pdf", plot = p, base_asp = 1/3, base_width = 9)
