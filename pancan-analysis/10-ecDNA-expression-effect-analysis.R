# Analyze in both gene and cytobands

library(gcap)
library(data.table)
library(IDConverter)
options(IDConverter.datapath = system.file("extdata", package = "IDConverter"))

setwd(file.path(PROJ_DIR, "pancan-analysis"))

proj_files = list.files("/data3/wsx_data/tcga_snp_array_gcap_result", pattern = "prediction_result", full.names = TRUE, all.files = TRUE)
fcna_list = purrr::map_df(proj_files, readRDS)
fcna_list
fcna = fcna_list[, .(sample, gene_id, total_cn, prob)]
rm(fcna_list); gc()
fcna[, is_circular := data.table::fifelse(prob > 0.5, 1L, 0L)] # 直接采用预测模型结果分组，更为宽松

fcna_gene = fcna[, .(N = sum(is_circular, na.rm = TRUE)), by = .(gene_id)]
fcna_gene = fcna_gene[N > 1]

fcna = fcna[gene_id %in% fcna_gene$gene_id]
fcna

fcna_ids = unique(fcna$sample)
fcna_ids_short = substr(fcna_ids, 1, 15)
sum(duplicated(fcna_ids_short))

# Remove repeated samples, as they are ambiguous
fcna_ids = fcna_ids[!fcna_ids_short %in% fcna_ids_short[duplicated(fcna_ids_short)]]
fcna = fcna[sample %in% fcna_ids]
fcna[, sample := substr(sample, 1, 15)]

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Prepare Gene Expression data
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
expr = data.table::fread("/data3/wsx_data/Xena/tcga_RSEM_gene_tpm.gz")
expr[, sample := substr(sample, 1, 15)]

sum(duplicated(expr$sample))

expr2 = melt(expr, id.vars = "sample")[
  !is.na(sample), .(value = mean(value, na.rm = TRUE)), by = .(sample, variable)] %>% 
  dcast(variable ~ sample)
colnames(expr2)[1] = "sample"
# log2(tpm+0.001)

expr2[1:5, 1:5]
ids = intersect(unique(fcna$sample), expr2$sample)
fcna = fcna[sample %in% ids]
expr2 = expr2[sample %in% ids]

fcna_gene = fcna[, .(N = sum(is_circular, na.rm = TRUE)), by = .(gene_id)]
fcna_gene = fcna_gene[N > 2]

fcna_gene[order(N)]
fcna = fcna[gene_id %in% fcna_gene$gene_id]
fcna

ids = fcna_gene$gene_id[fcna_gene$gene_id %in% colnames(expr2)]

expr = expr2[, c("sample", ids), with = FALSE]
expr = expr[sample %in% unique(fcna$sample)]
expr[1:5, 1:5]
expr_dt = melt(expr, id.vars = "sample", variable.name = "gene_id", value.name = "tpm")
expr_dt

mdt = merge(fcna, expr_dt, by = c("sample", "gene_id"))
mdt

saveRDS(mdt, file = "/data3/wsx_data/Xena/CN_Expr_fCNA_merged.rds")
rm(list = ls()); gc()

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Analysis
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
data = readRDS("/data3/wsx_data/Xena/CN_Expr_fCNA_merged.rds")

library(UCSCXenaShiny)

data2 = merge(data, 
              tcga_clinical %>%
                dplyr::select(sample, type) %>% 
                unique() %>% 
                data.table::as.data.table(), by = "sample", all.x = TRUE)
sum(is.na(data2$type))
data = data2; rm(data2); gc()
data[, type := factor(type)]
data = data[order(gene_id)]

genes = unique(data$gene_id)

gene = genes[1]
dt = data[gene_id %in% gene]
dt

# plot(tpm ~ total_cn, data =  dt[total_cn > 6])
# abline(coef = fit$coefficients)

analyzeLModel = function(dt) {
  # TODO: 加上circular/noncircular比值
  g = dt$gene_id[1]
  message("Processing gene ", g)
  fit = lm(tpm ~ total_cn + is_circular + type, data = dt)
  
  mresult = broom::tidy(fit)
  
  rv = data.table::data.table(
    gene_id = g,
    estimate_CN = mresult$estimate[mresult$term == "total_cn"],
    pvalue_CN = mresult$p.value[mresult$term == "total_cn"],
    estimate_Circular = mresult$estimate[mresult$term == "is_circular"],
    pvalue_Circular = mresult$p.value[mresult$term == "is_circular"],
    CN_mean_circular = mean(dt$total_cn[dt$is_circular == 1], na.rm = TRUE),
    Nobs = broom::glance(fit)$nobs,
    Ncircular = sum(dt$is_circular)
  )
  print(rv)
  rv
}

analyzeLModel(dt)

#data2 = data[order(gene_id)][1:20000]
result_allsamps = purrr::map_df(split(data, data$gene_id), analyzeLModel)
result_allsamps$p_adj_CN = p.adjust(result_allsamps$pvalue_CN, method = "BH")
result_allsamps$p_adj_Circular = p.adjust(result_allsamps$pvalue_Circular, method = "BH")
result_allsamps$symbol = IDConverter::convert_hm_genes(result_allsamps$gene_id)
result_allsamps$estimate_ratio = result_allsamps$estimate_Circular / result_allsamps$estimate_CN



saveRDS(result_allsamps, file = "data/Circular_Expr_effects.rds")

# Only plot biological sound data
data = result_allsamps[estimate_CN > 0 & estimate_Circular > 0][order(p_adj_Circular, -Ncircular)]
openxlsx::write.xlsx(
  data,
  file = "data/Circular_Expr_effects.xlsx"
)

library(ggplot2)
p = ggplot(subset(data, p_adj_Circular >= 0.05), aes(Ncircular, -log10(p_adj_Circular))) +
  geom_point() +
  geom_point(aes(color = estimate_Circular), data = subset(data, p_adj_Circular < 0.05), color = "steelblue") +
  ggrepel::geom_label_repel(
    aes(label = symbol), data = subset(
      data, (Ncircular > 110 & p_adj_Circular < 0.05) | -log10(p_adj_Circular) > 20),
    max.overlaps = 100
  ) +
  ggrepel::geom_label_repel(
    aes(label = symbol), data = subset(
      data, Ncircular > 150 & p_adj_Circular >= 0.05),
    max.overlaps = 100, color = "red", min.segment.length = 2
  ) +
  cowplot::theme_cowplot() +
  labs(x = "Frequency", y = "-log10(p.adjust)")
p
ggsave(filename = "plots/Circular_Expr_effective_geneplot.pdf",plot = p, width = 8, height = 7)


# MYC 不显著
# 是不是FGFR2的转录影响MYC转录导致的？
data3 = readRDS("/data3/wsx_data/Xena/CN_Expr_fCNA_merged.rds")
# MYC ENSG00000136997
# FGFR2 ENSG00000066468
data2 = data3[gene_id %in% c("ENSG00000136997", "ENSG00000066468")]
data2

plot(tpm ~ total_cn, data = data2[gene_id == "ENSG00000136997"])
plot(tpm ~ total_cn, data = data2[gene_id == "ENSG00000066468"])

dt = merge(data2[gene_id == "ENSG00000136997"][, .(sample, tpm, CN_MYC = total_cn, Circular_MYC = is_circular)],
           data2[gene_id == "ENSG00000066468"][, .(sample, Expr_FGFR2 = tpm, CN_FGFR2 = total_cn, Circular_FGFR2 = is_circular)],
           by = "sample")
dt

dt = merge(dt, 
           tcga_clinical %>%
             dplyr::select(sample, type) %>% 
             unique() %>% 
             data.table::as.data.table(), by = "sample", all.x = TRUE)
dt[, type := factor(type)]
dt

fit = lm(tpm ~ CN_MYC + Circular_MYC + Expr_FGFR2 + type, data = dt)
summary(fit)
p = forestmodel::forest_model(fit, covariates = c("CN_MYC", "Circular_MYC", "Expr_FGFR2"))
ggsave("plots/FRFR2_and_MYC_lm_forest.pdf", width = 8, height = 3)

mresult = broom::tidy(fit)

# 看下ERBB2
# ENSG00000141736
data4 = data3[gene_id %in% c("ENSG00000141736")]
data4

data4 = merge(data4, 
      tcga_clinical %>%
        dplyr::select(sample, type) %>% 
        unique() %>% 
        data.table::as.data.table(), by = "sample", all.x = TRUE)
data4[, type := factor(type)]


fit1 = lm(tpm ~ total_cn + is_circular + type, data = data4)
fit2 = lm(tpm ~ total_cn + type, data = data4)

dt = data.table::copy(data4)
dt[, `:=`(y_fit1 = as.numeric(predict(fit1, data4)), y_fit2 = as.numeric(predict(fit2, data4)))]
dt

p = ggplot(dt, aes(total_cn, tpm)) +
  geom_point() +
  geom_line(aes(y = y_fit1), color = "red") +
  geom_line(aes(y = y_fit2), color = "blue") +
  cowplot::theme_cowplot()
ggsave("plots/ERBB2.pdf", width = 5, height = 4)
