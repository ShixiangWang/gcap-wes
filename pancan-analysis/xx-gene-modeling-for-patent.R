library(dplyr)

data = readRDS("modeling/data/train_data_v3_raw.rds")
data
load("modeling/data/train_data_extra_info_v3.RData")

data = data[, .(sample, gene_id, total_cn = segVal, ecStatus = y)]

load("/data3/wsx_data/tcga.RData")

get_data = function(x) {
  data.table::melt(tcga[[x]][sample %in% data_gene_load$sample], id.vars = "sample")
}

expr = get_data("expression")
colnames(expr)[2:3] = c("gene_id", "expression")

mutation = get_data("mutation")
colnames(mutation)[2:3] = c("gene_id", "mutation")

methylation = get_data("methylation")
colnames(methylation)[2:3] = c("gene_id", "methylation")

fusion = rbind(tcga$fusion[, c(2, 3)] %>% setNames(c("sample", "gene_id")),
               tcga$fusion[, c(2, 4)] %>% setNames(c("sample", "gene_id")))
fusion = subset(fusion, sample %in% data_gene_load$sample)
fusion$fusion = 1L
fusion = unique(fusion) %>% data.table::as.data.table()

refs = readRDS("preprocessing-and-EDA/data/hg38_gene_info.rds")
refs[, gene_id := substr(gene_id, 1, 15)]
fusion[, gene_id := IDConverter::convert_custom(gene_id, "gene_name", "gene_id", refs)]
fusion = fusion[!is.na(gene_id)]

data2 = purrr::reduce(list(
  data,
  expr,
  mutation,
  methylation,
  fusion
), data.table::merge.data.table, by = c("sample", "gene_id"), all.x = TRUE)

rm(tcga); gc()

data2[, fusion := ifelse(is.na(fusion), 0L, 1L)]
data2

# Is it possible to construct a model from the gene quantification?
data3 = na.omit(data2)
data3[, cn := fcase(
  total_cn > 2, "AMP",
  total_cn == 2, "NORMAL",
  total_cn < 2, "DEL"
)]

data3$cn = factor(data3$cn)
data3$cn = relevel(data3$cn, "NORMAL")

lm_fit = glm(ecStatus ~ cn + expression + mutation + methylation + fusion, data = data3, family = binomial())
summary(lm_fit)

# baseline https://stats.stackexchange.com/questions/251175/what-is-baseline-in-precision-recall-curve

# Combined
pred = predict(lm_fit, type = "response")

get_auc = function(x, y) {
  gcap::get_auc(x, y)$auc.integral
}

# Save
elementary_predictors = data.frame(
  baseline = sum(data2$ecStatus) / nrow(data2),
  cn_DEL = get_auc(ifelse(data2$total_cn < 2, 1, 0), data2$ecStatus),
  cn_AMP = get_auc(ifelse(data2$total_cn > 2, 1, 0), data2$ecStatus),
  expression = get_auc(data2$expression[!is.na(data2$expression)], data2$ecStatus[!is.na(data2$expression)]),
  mutation = get_auc(data2$mutation[!is.na(data2$mutation)], data2$ecStatus[!is.na(data2$mutation)]),
  methylation = get_auc(data2$methylation[!is.na(data2$methylation)], data2$ecStatus[!is.na(data2$methylation)]),
  fusion = get_auc(data2$fusion, data2$ecStatus),
  combined = get_auc(pred, data3$ecStatus)
)

# Visualization
library(ggplot2)

data = elementary_predictors
colnames(data)[8] = "combined\nwith logistic\nregression"
data.table::setDT(data)
data = data.table::melt(data, measure.vars = colnames(data), variable.name = "model", value.name = "auPR")

p = ggplot(data, aes(x = model, y = auPR)) +
  geom_bar(stat = "identity", width = 0.6) +
  geom_text(aes(label = round(auPR, 3)), nudge_y = 0.003, color = "blue") +
  cowplot::theme_cowplot() +
  ggpubr::rotate_x_text(45) +
  labs(x = NULL, y = "PR-AUC")
p

ggsave("pancan-analysis/plots/auPR_of_elementary_predictors_4_patent.pdf", 
       p,
       width = 6, height = 4)
