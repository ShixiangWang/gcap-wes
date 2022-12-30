library(gcap)
library(data.table)

data = readRDS("modeling/data/train_data_v3_raw.rds")
colnames(data)[3] <- "total_cn"


result <- lapply(c("XGB11", "XGB32", "XGB56"),
                 function(x) gcap.runPrediction(data, model = x))
names(result) <- c("XGB11", "XGB32", "XGB56")


df <- cbind(as.data.table(result), data[, .(sample, y)])


# Overall performance -----------------------------------------------------
auc = list(
  auPRC = data.frame(
    XGB11 = get_auc(df$XGB11, df$y)$auc.integral,
    XGB32 = get_auc(df$XGB32, df$y)$auc.integral,
    XGB56 = get_auc(df$XGB56, df$y)$auc.integral
  ),
  auROC = data.frame(
    XGB11 = get_auc(df$XGB11, df$y, type = "roc")$auc,
    XGB32 = get_auc(df$XGB32, df$y, type = "roc")$auc,
    XGB56 = get_auc(df$XGB56, df$y, type = "roc")$auc
  )
)

auc

saveRDS(auc, file = "modeling/data/xgb_v3/overall_performance.rds")

# Cancer type importance --------------------------------------------------

df_all_cli <- readRDS("modeling/data/model_tcga_amplicon_related_clinical_data.rds")
df_all_cli %>%
  dplyr::ungroup() -> df_all_cli
df <- merge(df, df_all_cli %>%
              select(sample, type) %>%
              as.data.table(), by = "sample",
            all.x = TRUE)

## type特征的存在是否重要
pairs <- list(
  c("XGB32", "y"),
  c("XGB56", "y")
)

type_pair_results <- df[, lapply(pairs, function(x) {
  gcap::get_auc(.SD[[x[1]]], .SD[[x[2]]])$auc.integral
}), by = .(type)]

names(type_pair_results)[2:3] = c("XGB32", "XGB56")
type_pair_results = type_pair_results[!is.na(XGB56)]  # 3 type omited due to uncomputable AUC

type_pair_results[, `:=`(diff = XGB56 - XGB32)]
count_df <- unique(df[, .(type, sample)])[, .N, by = type]
type_pair_results <- merge(type_pair_results, count_df, by = "type", all.x = TRUE)

saveRDS(type_pair_results, file = "modeling/data/XGB32_56_type_importance.rds")

# Visualization
type_pair_results <- readRDS("modeling/data/XGB32_56_type_importance.rds")
summary(type_pair_results$diff)

library(ggplot2)
library(forcats)
library(lemon)

type_pair_results = type_pair_results[N >= 10]
dd = melt(type_pair_results, measure.vars = c("XGB32", "XGB56"))
#dd = dd[!type %in% "KIRC"]  # What's the reason causing bad performance on KIRC?
dd[, type := paste0(type, " (N=", N, ")")]
dd$type = fct_reorder(dd$type, dd$diff)

ggplot(data = dd, 
       mapping = aes(x = ifelse(variable == "XGB56", 
                                yes = value, no = -value), 
                     y = type, fill = variable)) +
  geom_col() +
  geom_text(aes(label = round(value, 2)), 
            data = dd[variable == "XGB32"],
            hjust = -1) +
  geom_text(aes(label = round(value, 2)), 
            data = dd[variable == "XGB56"],
            hjust = 1) +
  geom_text(aes(label = round(diff, 3)),
            data = dd[variable == "XGB56"], hjust = -0.2,
            color = "red") +
  scale_x_symmetric(labels = abs) +
  labs(x = "auPRC", y = NULL, fill = "Model") +
  cowplot::theme_cowplot() +
  theme(legend.position = "top") -> p
p

ggsave("modeling/plots/auPRC_comparison_across_types_N10.pdf", plot = p,
       width = 8, height = 6)

# Check the sample distribution in train and test sets
dd

source("lib/get_data.R")
dfs = get_model_data(only_samples = TRUE)
dfs2 = data.frame(
  set = c(rep("train", length(dfs$train)), rep("test", length(dfs$test))),
  sample = c(dfs$train, dfs$test)
)
dfs2 = dplyr::left_join(dfs2, df_all_cli, by = "sample")
dfs2

table(dfs2$set)
table(dfs2$set, dfs2$type)
# > table(dfs2$set, dfs2$type)
# 
#       BLCA BRCA CESC COAD DLBC ESCA GBM HNSC KICH KIRC KIRP LGG LIHC LUAD LUSC OV PRAD READ SARC SKCM STAD THCA UCEC UVM
# test     4    3    2    2    1    6   6    8    1    0    1   5    4    5    1  3    2    0    6    9    4    0    4   0
# train   29   26   18    5    0   26  20   35    1    2    1  10    4   20   14  5    7    1   11   19   27    5   20   3


