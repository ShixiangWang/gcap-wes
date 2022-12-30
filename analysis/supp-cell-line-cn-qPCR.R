setwd("manuscript/")
# 以CCLE标记为正常CN=2的作为参考基因

library(tidyverse)

data = read_csv("data/gene_copy_per_ng.csv")
data
unique(data$Gene)


# 以PC3的FGFR2作为参考
data = data %>% 
  mutate(CN = 2 * CN_mean / CN_mean[Gene == "FGFR2" & Sample == "PC3"])

data
ggplot(data %>% filter(Gene %in% c("MYC", "FGFR2")),
       aes(x = Sample, y = CN, fill = Gene)) +
  geom_bar(stat = "identity", position = position_dodge2(), width = 0.5) +
  cowplot::theme_cowplot() + theme(legend.position = "top") +
  labs(x = NULL, y = "Estimated copy number by qPCR") +
  geom_hline(yintercept = 2, linetype = 2) +
  coord_flip() -> p
p
ggsave("plots/SNU16_PC3_CN_estimated_by_qPCR.pdf", width = 5, height = 3)

# Compare CN estimated by qPCR and WES and WGS
data2 = fread("data/cell_line_gcap.tsv")
data2 = data2[gene_name %in% data$Gene & sample %in% c(
  "SNU16-1_WES", "pc3-1_WES"#, "pc3_WGS", "snu16_WGS",
)][, list(sample, gene_name, total_cn)]
data2

data2[, sample := ifelse(sample == "pc3-1_WES", "PC3", "SNU16")]
data2
data2$assay = "WES"

data3 = rbind(data2,
              data %>% 
                rename(gene_name = Gene,
                       sample = Sample,
                       total_cn = CN) %>% 
                mutate(assay = "qPCR") %>% 
                select(sample, gene_name, total_cn, assay) %>% 
                as.data.table())

data3

data3_wide = dcast(data3, sample + gene_name ~ assay, value.var = "total_cn")

p = ggplot(data3_wide, aes(x = WES, y = qPCR)) +
  geom_point() +
  geom_smooth(method = "lm") +
  ggpubr::stat_cor(method = "pearson", label.x.npc = 0.05, label.y.npc = 1) +
  facet_wrap(~sample) +
  scale_x_log10() +
  scale_y_log10() +
  ggrepel::geom_text_repel(aes(label = gene_name)) +
  #ggthemes::theme_base() +
  cowplot::theme_cowplot() +
  labs(x = "Copy number estimated by WES", y = "Copy number estimated by qPCR") +
  coord_equal()
p
ggsave("plots/SNU16_PC3_CN_estimated_by_qPCR_vs_WES.pdf", plot = p, width = 7, height = 4)
