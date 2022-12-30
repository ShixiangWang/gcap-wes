# Analyze in both gene and cytobands
PROJ_DIR = "~/gcap-analysis/manuscript/"
setwd(PROJ_DIR)

library(gcap)
library(data.table)
library(IDConverter)
options(IDConverter.datapath = system.file("extdata", package = "IDConverter"))

TCGA = readRDS("data/TCGA_SNP.rds")
TCGA$sample_summary = readRDS("data/tcga_snp_data.rds")

TCGA$getGeneSummary()

fcna = TCGA$data
# 初步筛选
fcna_gene = TCGA$getGeneSummary()[circular > 10]

# fcna = fcna[gene_id %in% fcna_gene$gene_id]
# fcna
# fcna_ids = unique(fcna$sample)

fcna[, gene_name := IDConverter::convert_hm_genes(gene_id, genome_build = "hg19")]
fcna_gene[, gene_name := IDConverter::convert_hm_genes(gene_id, genome_build = "hg19")]

fcna_gene[gene_name == "MYC"]
fcna_gene[gene_name == "EGFR"]
fcna_gene[gene_name == "ERBB2"]
fcna[gene_name == "MYC"][order(-total_cn)]

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Prepare Gene Expression data
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
expr = data.table::fread("/data3/wsx_data/Xena/tcga_RSEM_gene_tpm")
expr[, sample := substr(sample, 1, 15)]
olap_sample = colnames(expr)[colnames(expr) %in% TCGA$sample_summary$sample]
expr = expr[sample %in% unique(fcna$gene_id), c("sample", olap_sample), with = FALSE]
expr[1:5, 1:5]

sum(duplicated(expr$sample))

expr2 = melt(expr[sample %in% fcna_gene$gene_id], id.vars = "sample")[
  !is.na(sample), .(value = mean(value, na.rm = TRUE)), by = .(sample, variable)] %>% 
  dcast(variable ~ sample)
colnames(expr2)[1] = "sample"
expr2[1:5, 1:5]
# log2(tpm+0.001)

expr_dt = melt(expr2, id.vars = "sample", variable.name = "gene_id", value.name = "tpm")
expr_dt

proj_files = list.files("/data3/wsx_data/tcga_snp_array_gcap_result2/", 
                        pattern = "prediction_result", full.names = TRUE, all.files = TRUE)
cn_info = purrr::map_df(proj_files, function(x) {
  rv = readRDS(x)
  rv = rv[gene_id %in% fcna_gene$gene_id, 
          c("sample", "band", "gene_id", "total_cn", "prob", "gene_class"), with = FALSE]
  print(head(rv))
  rv
})
cn_info = cn_info[gene_id %in% unique(expr_dt$gene_id)]
cn_info = cn_info[sample %in% unique(expr_dt$sample)]

length(unique(expr_dt$gene_id))
length(unique(cn_info$gene_id))

length(unique(expr_dt$sample))
length(unique(cn_info$sample))

mdt = merge(cn_info, expr_dt, by = c("sample", "gene_id"))
mdt

saveRDS(mdt, file = "/data3/wsx_data/Xena/CN_Expr_fCNA_merged.rds")
rm(list = ls()); gc()

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Analysis
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
data = readRDS("/data3/wsx_data/Xena/CN_Expr_fCNA_merged.rds")

library(UCSCXenaShiny)

data2 = merge(data, 
              readRDS("data/tcga_snp_data.rds")[, list(sample, purity, type)], by = "sample")

data = data2; rm(data2); gc()
data[, type := factor(type)]
data = data[order(gene_id)]
data[, gene_name := IDConverter::convert_hm_genes(gene_id, genome_build = "hg19")]

genes = unique(data$gene_id)

analyzeLModel = function(dt) {
  g = dt$gene_id[1]
  #dt = dt[type %in% dt[, .N, by = list(gene_class, type)][gene_class == "circular"][N >= 3]$type]
  message("Processing gene ", g, " with name ", dt$gene_name[1])
  fit = lm(tpm ~ total_cn + is_circular + type + purity, data = dt)
  
  mresult = broom::tidy(fit)
  
  rv = data.table::data.table(
    gene_id = g,
    gene_name = dt$gene_name[1],
    pvalue_CN = mresult$p.value[mresult$term == "total_cn"],
    pvalue_Circular = mresult$p.value[mresult$term == "is_circular"],
    CN_mean_circular = mean(dt$total_cn[dt$is_circular == 1], na.rm = TRUE),
    CN_mean_nocircular = mean(dt$total_cn[dt$is_circular != 1], na.rm = TRUE),
    Nobs = broom::glance(fit)$nobs,
    Ncircular = sum(dt$is_circular)
  )
  print(rv)
  rv
}

data[, is_circular := ifelse(gene_class == "circular", 1L, 0L)]
data[, gene_class := gcaputils::set_default_factor(gene_class)]

ggstatsplot::ggbetweenstats(
  data[gene_name == "MYC"],
  x = gene_class,
  y = tpm, bf.message = FALSE
)

gene = genes[1]
dt = data[gene_id %in% gene]
dt
analyzeLModel(dt)

# Only focus on Oncogenes
#intogen = data.table::fread("~/gcap-analysis/manuscript/data/Compendium_Cancer_Genes.tsv")
#oncogenes = unique(intogen[ROLE == "Act"]$SYMBOL)
# data2 = data[gene_name %in% oncogenes]
# sum(unique(data$gene_name) %in% oncogenes)
# # 79

data2 = data[gene_name %in% gcap::oncogenes$OncogeneName]
sum(unique(data$gene_name) %in% gcap::oncogenes$OncogeneName)
# 236 最上面筛过一遍了

result_allsamps = purrr::map_df(split(data2, data2$gene_id), analyzeLModel)
result_allsamps$p_adj_CN = p.adjust(result_allsamps$pvalue_CN, method = "BH")
result_allsamps$p_adj_Circular = p.adjust(result_allsamps$pvalue_Circular, method = "BH")

TCGA = readRDS("data/TCGA_SNP.rds")
band_data = unique(TCGA$data[, .(band, gene_id)])
result = merge(result_allsamps, band_data, by = "gene_id", all.x = TRUE)
sum(is.na(result$band))

band_ord = result[, .(m = min(p_adj_Circular, na.rm = TRUE)), by = band][order(m)]
result[, band := factor(band, band_ord$band)]
result = result[order(p_adj_Circular)]
data.table::setcolorder(result, "band")

saveRDS(result, file = "data/Circular_Expr_effects2.rds")
if (file.exists("data/Circular_Expr_effects2.xlsx")) file.remove("data/Circular_Expr_effects2.xlsx")
openxlsx::write.xlsx(
  result[, list(band, gene_name, Nobs, p = pvalue_Circular, p_adj = p_adj_Circular)],
  file = "data/Circular_Expr_effects2.xlsx"
)

# # 基因数量太多了，我们用cytoband展示
# data2 = result[, .(p_adj_Circular = min(p_adj_Circular, na.rm = TRUE),
#                    N = sum(Ncircular, na.rm = TRUE),
#                    Size = sum(p_adj_Circular < 0.05, na.rm = TRUE)), by = band]


library(ggplot2)
p = ggplot(result, aes(CN_mean_nocircular, CN_mean_circular - CN_mean_nocircular)) +
  geom_point(aes(color = p_adj_Circular < 0.05)) +
  ggrepel::geom_label_repel(
    aes(label = gene_name), data = subset(
      result, p_adj_Circular < 0.05),
    max.overlaps = 100, size = 2.5
  ) +
  cowplot::theme_cowplot() +
  labs(y = "Copy number mean on ecDNA", x = "Copy number mean on chromosome",
       color = NULL) + theme(legend.position = "none") +
  scale_color_manual(values = c("grey", "red"))
p
ggsave(filename = "plots/Oncogene_Circular_Expr_effective_geneplot2.pdf",plot = p, width = 8, height = 7)

result2 = copy(result)
result2[, band2 := sub("(\\..+)", "", band)]
result3 = result2 %>%
  dplyr::group_by(band2) %>%
  dplyr::top_n(-p_adj_Circular, n = 1)

p = ggplot(result3, aes(CN_mean_nocircular, CN_mean_circular - CN_mean_nocircular)) +
  geom_point(aes(color = p_adj_Circular < 0.05)) +
  ggrepel::geom_label_repel(
    aes(label = band2), data = subset(
      result3, p_adj_Circular < 0.05),
    max.overlaps = 100, size = 2.5
  ) +
  cowplot::theme_cowplot() +
  labs(y = "Copy number mean on ecDNA", x = "Copy number mean on chromosome",
       color = NULL) + theme(legend.position = "none") +
  scale_color_manual(values = c("grey", "red"))
p
ggsave(filename = "plots/Oncogene_Circular_Expr_effective_band_most_sig_geneplot.pdf",plot = p, width = 8, height = 7)


plot_genes = result[p_adj_Circular < 0.05]$gene_name
library(ggpubr)
p = ggboxplot(data[gene_name %in% plot_genes], x = "gene_class", y = "tpm", fill = "gene_class",
              facet.by = "gene_name", scales = "free_y", xlab = FALSE, ylab = "Gene expression (log2 based TPM)") +
  rotate_x_text(45) + theme(legend.position = "none")
ggsave(filename = "plots/Oncogene_fCNA_class_expression2.pdf",plot = p, width = 18, height = 15)

# Fold change LM like NG
dt = copy(data)
dt[, tpm := 2^tpm - 0.001]
# plus 1 to avoid 0
da_up_nofocal = dt[gene_class == "nofocal", list(ref = mean(tpm[tpm > quantile(tpm, 0.75)], na.rm = TRUE) + 1),
                      by = list(gene_name)]
da_up25_circular = dt[gene_class == "circular"][, .SD[tpm > quantile(tpm, 0.75)], by = "gene_name"]
da_up25_noncircular = dt[gene_class == "noncircular"][, .SD[tpm > quantile(tpm, 0.75)], by = "gene_name"]
da = merge(rbind(da_up25_circular, da_up25_noncircular),
           da_up_nofocal, by = "gene_name")
da[, foldchange := (tpm + 1) / ref]

p = ggplot(da, aes(x = total_cn, y = foldchange, color = gene_class)) +
  geom_smooth(method = "lm", fullrange = TRUE) +
  labs(x = "Copy number", y = "TPM upper quantile", color = "Classification") +
  cowplot::theme_cowplot() +
  scale_color_manual(values = c("#0066CC", "#CC0033"))
ggsave("plots/TCGA_TPM_foldchange_vs_CN2.pdf", plot = p, width = 5, height = 4)

dt = copy(data2)
dt[, tpm := 2^tpm - 0.001]
# plus 1 to avoid 0
da_up_nofocal = dt[gene_class == "nofocal", list(ref = mean(tpm[tpm > quantile(tpm, 0.75)], na.rm = TRUE) + 1),
                   by = list(gene_name)]
da_up25_circular = dt[gene_class == "circular"][, .SD[tpm > quantile(tpm, 0.75)], by = "gene_name"]
da_up25_noncircular = dt[gene_class == "noncircular"][, .SD[tpm > quantile(tpm, 0.75)], by = "gene_name"]
da = merge(rbind(da_up25_circular, da_up25_noncircular),
           da_up_nofocal, by = "gene_name")
da[, foldchange := (tpm + 1) / ref]

p = ggplot(da, aes(x = total_cn, y = foldchange, color = gene_class)) +
  geom_smooth(method = "lm", fullrange = TRUE) +
  labs(x = "Copy number", y = "TPM upper quantile", color = "Classification") +
  cowplot::theme_cowplot() +
  scale_color_manual(values = c("#0066CC", "#CC0033"))
ggsave("plots/TCGA_TPM_foldchange_vs_CN_only_known_oncogenes2.pdf", plot = p, width = 5, height = 4)

# all available genes

data_gene = data[gene_name %in% result$gene_name]
data_gene[, class := ifelse(is_circular, "circular", "nofocal/noncircular")]
data_gene[, class := factor(class, c("nofocal/noncircular", "circular"))]

library(gghalves)

p_stat = ggpubr::compare_means(tpm ~ class, data = data_gene, method = "t.test", group.by = "gene_name")

kept_genes = p_stat$gene_name[p_stat$p.adj < 0.05]
# too many (176), plot top50
p_stat2 = p_stat %>% dplyr::arrange(p.adj) %>% head(50)

p = ggplot2::ggplot(
  data_gene[gene_name %in% p_stat2$gene_name] %>%
    dplyr::mutate(gene_name = factor(gene_name, p_stat2$gene_name))
) +
  gghalves::geom_half_violin(aes(
    x = gene_name, y = tpm, fill = class, color = class,
    split = class), position = "identity") +
  ggplot2::xlab(NULL) + ylab("Gene expression (log2 based TPM)") +
  ggplot2::scale_fill_manual(values = c("grey", "red")) +
  ggplot2::scale_color_manual(values = c("grey", "red")) +
  cowplot::theme_cowplot() +
  rotate_x_text(45) +
  theme(axis.text.x = element_text(size = 8), legend.position = "top")

p
ggsave(filename = "plots/Oncogene_fCNA_class_expression_all_sigs_half_violin_top50.pdf",plot = p, width = 12, height = 4)
