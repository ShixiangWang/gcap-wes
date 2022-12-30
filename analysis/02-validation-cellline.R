PROJ_DIR = "~/gcap-analysis/manuscript/"
setwd(PROJ_DIR)

data = fread("data/cell_line_gcap.tsv")

# similarity
data_cn = na.omit(dcast(data[gene_name != ""], gene_name ~ sample,
                        value.var = "total_cn", fun.aggregate = function(x) x[1]))
p = ggcorrplot::ggcorrplot(cor(data_cn[, -1]), hc.order = TRUE, lab = TRUE)
p
ggsave("plots/cell_line_CN_profile_corrplot.pdf", plot = p, width = 7, height = 6)
nrow(data_cn)
# 18805

# Distribution
devtools::load_all("~/proj/gcaputils/")
p = gcap.plotDistribution(data[gene_class != "nofocal", 
                           list(sample = paste0(gene_name, sample),
                                class = gene_class, by = sample)],
                      fill = FALSE, 
                      palette = c("#0066CC", "#CC0033"))
p
ggsave("plots/cell_line_amplicon_dist.pdf", plot = p, width = 5, height = 4)


# WGS overlap
data_wgs = data[endsWith(sample, "WGS")]
data_wgs[, sample := stringr::str_remove(sample, "_WGS")]

# 虽然aa的结果有有一个文件report了Gene list，但和最后的summary结果不一致，所以只分析看能否找到summary
# 提到的基因
data_aa = readxl::read_excel("data/cellline_wgs_result.xlsx", sheet = 1)
data_aa = as.data.table(data_aa)
data_aa = data_aa[amplicon_decomposition_class != "No amp/Invalid", 
                  list(sample = id, AmpliconID, `ecDNA+`, amplicon_decomposition_class, Intervals, OncogenesAmplified)]
data_aa[, gene_class := ifelse(`ecDNA+` == "Positive", "circular", "noncircular")]

data_aa = data_aa %>% 
  tidyr::separate_rows(Intervals, sep = ",") %>% 
  dplyr::select(sample, Intervals, gene_class) %>% 
  tidyr::separate(Intervals, into = c("chr", "start", "end"), sep = "[:-]") %>%
  dplyr::mutate(start = as.integer(start), end = as.integer(end)) %>%
  data.table::as.data.table()

data.table::setcolorder(data_aa, c("chr", "start", "end", "gene_class", "sample"))

ref_genes <- readRDS("../preprocessing-and-EDA/data/hg38_gene_info.rds")
ref_genes <- ref_genes[chrom %in% paste0("chr", 1:22) & gene_type == "protein_coding"]
ref_genes$strand <- NULL
ref_genes[, gene_id := gsub("(\\..+)", "", gene_id)]
ref_genes$gene_type <- NULL

# Use gene_id is more robust
# > length(unique(ref_genes$gene_id))
# > length(unique(ref_genes$gene_name))
ref_genes$gene_name <- NULL
colnames(ref_genes)[1] = "chr"

data_wgs2 = merge(data_wgs, ref_genes, by = "gene_id")
data.table::setcolorder(data_wgs2, c("chr", "start", "end", "sample", "gene_id"))
data_wgs2 = data_wgs2[gene_class != "nofocal"][, list(chr, start, end, sample, band, gene_id, total_cn, prob, gene_class)]

data.table::setkey(data_wgs2, NULL)
regions <- modules::use("../lib/regions.R")
amp_dt = list()
for (i in unique(data_wgs2$sample)) {
  #undebug(regions$overlaps)
  amp_dt[[i]] = regions$overlaps(data_aa[sample == i, list(chr, start, end, aa_class = gene_class)],
                                 data_wgs2[sample == i])
}
amp_dt = rbindlist(amp_dt)

summary(amp_dt$intersect_ratio)
# # Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.4515  1.0000  1.0000  0.9970  1.0000  1.0000 

# The calculation is limited to overlap genes

amp_overlap = rbind(
  data.table(
    type = "circular", 
    score = c("precision", "recall"),
    pct = c(
      nrow(amp_dt[aa_class == "circular" & gene_class == "circular"]) / nrow(amp_dt[gene_class == "circular"]),
      nrow(amp_dt[aa_class == "circular" & gene_class == "circular"]) / nrow(amp_dt[aa_class == "circular"])
    )),
  data.table(
    type = "noncircular", 
    score = c("precision", "recall"),
    pct = c(
      nrow(amp_dt[aa_class == "noncircular" & gene_class == "noncircular"]) / nrow(amp_dt[gene_class == "noncircular"]),
      nrow(amp_dt[aa_class == "noncircular" & gene_class == "noncircular"]) / nrow(amp_dt[aa_class == "noncircular"])
    ))
)
amp_overlap

amp_overlap2 = list()
for (i in unique(amp_dt$sample)) {
  dt1 = amp_dt[sample == i]
  dt2 = data_wgs2[sample == i]
  dt3 = rbind(
    data.table(
      type = "circular", 
      score = c("precision", "recall"),
      pct = c(
        nrow(dt1[aa_class == "circular" & gene_class == "circular"]) / nrow(dt1[gene_class == "circular"]),
        nrow(dt1[aa_class == "circular" & gene_class == "circular"]) / nrow(dt1[aa_class == "circular"])
      )),
    data.table(
      type = "noncircular", 
      score = c("precision", "recall"),
      pct = c(
        nrow(dt1[aa_class == "noncircular" & gene_class == "noncircular"]) / nrow(dt1[gene_class == "noncircular"]),
        nrow(dt1[aa_class == "noncircular" & gene_class == "noncircular"]) / nrow(dt1[aa_class == "noncircular"])
      ))
  )
  amp_overlap2[[i]] = dt3
}
amp_overlap2 = rbindlist(amp_overlap2, idcol = "sample")
amp_overlap2 = amp_overlap2[!is.na(pct) & !(sample %in% c("HCC827", "HK301", "CAKI2", "T47D") & type == "circular")]

amp_overlap2
library(ggbeeswarm)
p = ggplot(amp_overlap2 %>% 
             dplyr::mutate(type = factor(type, c("circular", "noncircular"))), 
           aes(x = type, y = pct)) +
  geom_boxplot(width = 0.5, outlier.shape = NA) +
  geom_quasirandom(aes(color = sample), width = 0.2) +
  #geom_point(aes(x = type, y = pct), data = amp_overlap, color = "red", shape = 17, size = 4, alpha = 0.6) +
  facet_wrap(~score) +
  cowplot::theme_cowplot() +
  labs(x = NULL, y = NULL) + 
  scale_color_brewer(palette = "Set1")
p

ggsave("plots/cell_line_amplicon_ratio.pdf", plot = p, width = 6, height = 4)
save(amp_overlap, amp_overlap2, file = "data/cellline_amp_overlap.RData")

# PC3/SNU16 WGS and WES comparison
data2 = data[startsWith(sample, "pc3") | startsWith(sample, "SNU16") | startsWith(sample, "snu16")]
data2
data2 = data2 %>% tidyr::separate("sample", c("sample", "assay"), sep = "_")
data2[, sample := sub("-1", "", sample)]
data2[, sample := ifelse(sample == "snu16", "SNU16", sample)]
data2[sample == "pc3"]$assay %>% table()
data2[sample == "SNU16"]$assay %>% table()

dd_class = dcast(data2[!is.na(gene_name) & gene_name != ""], 
                 sample + gene_name ~ assay, value.var = "gene_class", fun.aggregate = function(x) x[1])
dd_class = dd_class[!is.na(WES) & !is.na(WGS)]
dd_cn = dcast(data2[!is.na(gene_name) & gene_name != ""], sample + gene_name ~ assay, value.var = "total_cn", fun.aggregate = mean)
dd_cn = dd_cn[!is.na(WES) & !is.na(WGS)]
dd = merge(dd_class, dd_cn)
dd[, type := ifelse(WES.x == WGS.x, "consistent", "inconsistent")]
dd = dd[!(WES.x == "nofocal" & WGS.x == "nofocal") & !is.na(WES.x) & !is.na(WGS.x)]
dd

ggplot(dd, aes(log2(WES.y), log2(WGS.y))) +
  geom_jitter(aes(color = type), alpha = 0.5, width = 0.5) +
  ggrepel::geom_text_repel(aes(label = gene_name), data = subset(dd, type == "consistent"), max.overlaps = 20,
                           size = 2) +
  facet_wrap(~sample, scales = "free_y") + 
  cowplot::theme_cowplot() +
  labs(x = "CN by WES (log2 based)", y = "CN by WGS (log2 based)") +
  theme(legend.position = "top") -> p

ggsave("plots/cell_line_pc3_and_snu16_amplicon_gene_comparison.pdf", p, width = 7,height = 4)


dd2 = merge(dd, unique(data2[, list(band, gene_name)]), by = "gene_name")

r1 = dd2[, .N, by = list(sample, type)]

r2 = dd2[type == "inconsistent"][, .N, by = list(sample, band)][order(-N)]
r1
r2

openxlsx::write.xlsx(
  list(consistence = r1,
       inconsistence_dist = r2),
  file = "data/pc3_snu16_consistence.xlsx"
)
