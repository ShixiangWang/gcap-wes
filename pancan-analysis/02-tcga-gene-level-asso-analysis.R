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
# > nrow(data3)
# [1] 4752969
lm_fit = glm(ecStatus ~ total_cn + expression + mutation + methylation + fusion, data = data3, family = binomial())
summary(lm_fit)

# Coefficients:
#   Estimate Std. Error  z value Pr(>|z|)    
# (Intercept) -10.365674   0.035309 -293.569  < 2e-16 ***
#   total_cn      0.943723   0.003612  261.294  < 2e-16 ***
#   expression    0.016831   0.002347    7.171 7.46e-13 ***
#   mutation      0.350961   0.092268    3.804 0.000143 ***
#   methylation  -0.343465   0.054424   -6.311 2.77e-10 ***
#   fusion        2.849211   0.145366   19.600  < 2e-16 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

# Use features can be obtained from RNAseq
# data4 = na.omit(data2[, .(expression, fusion, ecStatus)])
# lm_fit2 = glm(ecStatus ~ expression + fusion, data = data4, family = binomial())
# summary(lm_fit2)

library(ggeffects)
library(parameters)

mm = model_parameters(lm_fit, ci_method="wald", exponentiate = FALSE)
mm
# Parameter   | Log-Odds |       SE |           95% CI |       z |      p
# -----------------------------------------------------------------------
#   (Intercept) |   -10.37 |     0.04 | [-10.43, -10.30] | -293.57 | < .001
# total cn    |     0.94 | 3.61e-03 | [  0.94,   0.95] |  261.29 | < .001
# expression  |     0.02 | 2.35e-03 | [  0.01,   0.02] |    7.17 | < .001
# mutation    |     0.35 |     0.09 | [  0.17,   0.53] |    3.80 | < .001
# methylation |    -0.34 |     0.05 | [ -0.45,  -0.24] |   -6.31 | < .001
# fusion      |     2.85 |     0.15 | [  2.56,   3.13] |   19.60 | < .001
plot(mm)

library(ggplot2)
gg = plot(mm)
ggsave("pancan-analysis/plots/logistic_association_between_molecular_profile_and_ecStatus.pdf", gg,
       width = 6, height = 4)

table(ecStatus = data2$ecStatus, mutation = data2$mutation)
#     0       1
# 0 6071227   49877
# 1   21368     334
# ecStatus为1且突变的其实很少

fisher.test(table(ecStatus = data2$ecStatus, mutation = data2$mutation))
# Fisher's Exact Test for Count Data
# 
# data:  table(ecStatus = data2$ecStatus, mutation = data2$mutation)
# p-value < 2.2e-16
# alternative hypothesis: true odds ratio is not equal to 1
# 95 percent confidence interval:
#  1.702070 2.120709
# sample estimates:
# odds ratio 
#   1.902662 
# 接近2倍的富集

table(ecStatus = data2$ecStatus, fusion = data2$fusion)
# fisher.test(table(ecStatus = data2$ecStatus, fusion = data2$fusion))
# 
# Fisher's Exact Test for Count Data
# 
# data:  table(ecStatus = data2$ecStatus, fusion = data2$fusion)
# p-value < 2.2e-16
# alternative hypothesis: true odds ratio is not equal to 1
# 95 percent confidence interval:
#  46.72503 58.77413
# sample estimates:
# odds ratio 
#   52.42782 
# 52倍的富集

# baseline https://stats.stackexchange.com/questions/251175/what-is-baseline-in-precision-recall-curve

# Combined
pred = predict(lm_fit, type = "response")

get_auc = function(x, y) {
  gcap::get_auc(x, y)$auc.integral
}

# Save
elementary_predictors = data.frame(
  baseline = sum(data2$ecStatus) / nrow(data2),
  total_cn = get_auc(data2$total_cn, data2$ecStatus),
  expression = get_auc(data2$expression[!is.na(data2$expression)], data2$ecStatus[!is.na(data2$expression)]),
  mutation = get_auc(data2$mutation[!is.na(data2$mutation)], data2$ecStatus[!is.na(data2$mutation)]),
  methylation = get_auc(data2$methylation[!is.na(data2$methylation)], data2$ecStatus[!is.na(data2$methylation)]),
  fusion = get_auc(data2$fusion, data2$ecStatus),
  combined = get_auc(pred, data3$ecStatus)
)

elementary_predictors_n = data.frame(
  baseline = nrow(data2),
  total_cn = length(data2$total_cn),
  expression = length(data2$expression[!is.na(data2$expression)]),
  mutation = length(data2$mutation[!is.na(data2$mutation)]),
  methylation = length(data2$methylation[!is.na(data2$methylation)]),
  fusion = sum(!is.na(data2$fusion)),
  combined = nrow(data3)
)

save(elementary_predictors, elementary_predictors_n, file = "pancan-analysis/data/estimates_predictors.RData")

# Visualization
library(ggplot2)

load(file = "pancan-analysis/data/estimates_predictors.RData")
data = elementary_predictors
colnames(data)[2] = "copy number"
colnames(data)[7] = "combined\nwith logistic\nregression"
data.table::setDT(data)
data = data.table::melt(data, measure.vars = colnames(data), variable.name = "model", value.name = "auPR")

p = ggplot(data, aes(x = model, y = auPR)) +
  geom_bar(stat = "identity", width = 0.6) +
  geom_text(aes(label = round(auPR, 3)), nudge_y = 0.03, color = "blue") +
  cowplot::theme_cowplot() +
  ggpubr::rotate_x_text(45) +
  labs(x = NULL)
p

ggsave("pancan-analysis/plots/auPR_of_elementary_predictors.pdf", 
       p,
       width = 6, height = 4)

p = ggplot(data[!grepl("combined", model)], aes(x = model, y = auPR)) +
  geom_bar(stat = "identity", width = 0.6) +
  geom_text(aes(label = round(auPR, 3)), nudge_y = 0.03, color = "blue") +
  cowplot::theme_cowplot() +
  ggpubr::rotate_x_text(45) +
  labs(x = NULL, y = "auPRC")
p

ggsave("pancan-analysis/plots/auPR_of_elementary_predictors2.pdf", 
       p,
       width = 6, height = 4)

