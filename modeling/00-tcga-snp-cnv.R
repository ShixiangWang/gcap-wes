setwd(file.path(PROJ_DIR, "modeling"))

library(dplyr)

data_snp = readRDS("../preprocessing-and-EDA/data/tcga_ascat.rds")

ploidy = unique(data_snp[, .(sample, ploidy)])

hist(ploidy$ploidy, breaks = 100)
# many samples' ploidy over 3

ploidy_comps = gcaputils::gcap.extractComponents(ploidy$ploidy, threshold = 0.3)
ploidy_comps$params

# > ploidy_comps$params
# mean        sd    n
# 1: 2.003470 0.1195671 4867
# 2: 3.518504 0.7720292 4832

# Only keep sample with ploidy within 1SD
# ploidy = ploidy[ploidy < 2.003470 + 0.1195671 & ploidy > 2.003470 - 0.1195671]
# ploidy
#saveRDS(ploidy, file = "data/selected_samples.rds")

data_gene = gcap:::collapse_to_genes(data_snp[, .(chr, start, end, sample, total_cn)], drop = FALSE)
data_gene = data_gene[intersect_ratio == 1]


# Drop out outliers to obtain a robust estimation of mean and sd for gene copy number profile
gene_cn = data_gene[, .(
  somatic_cn_mean = mean(total_cn[!rstatix::is_extreme(total_cn)], na.rm = TRUE),
  somatic_cn_sd = sd(total_cn[!rstatix::is_extreme(total_cn)], na.rm = TRUE)
), by = .(gene_id)]

summary(gene_cn$somatic_cn_mean)
summary(gene_cn$somatic_cn_sd)

hist(gene_cn$somatic_cn_mean, breaks = 100)
hist(gene_cn$somatic_cn_sd, breaks = 100)

options(IDConverter.datapath = system.file("extdata", package = "IDConverter"))
gene_cn$gene_name = IDConverter::convert_hm_genes(gene_cn$gene_id)
gene_cn

gene_cn = gene_cn[order(somatic_cn_mean)]

saveRDS(gene_cn, file = "data/somatic_gene_cn.rds")
