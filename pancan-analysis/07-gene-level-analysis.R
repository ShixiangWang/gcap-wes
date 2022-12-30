setwd(file.path(PROJ_DIR, "pancan-analysis"))

TCGA = readRDS("data/TCGA_SNP.rds")
# TCGA Gene analysis -----------------------------------------------------------

TCGA$convertGeneID(genome_build = "hg38")
gene_info = TCGA$gene_summary[!is.na(gene_id)][order(circular, Total, decreasing = TRUE)]
sum(gene_info$circular > 1)

#gcap.plotProfile(TCGA, top_n = 100, show_column_names = FALSE)
#gcap.plotDistribution(TCGA, x = gene_info[circular > 2]$gene_id)

gcap.plotCircos(TCGA, genome_build = "hg38")

genes_summary = TCGA$data[
  , .(cn = mean(total_cn[amplicon_type %in% c("circular", "possibly_circular")], na.rm = TRUE),
      prob = mean(prob[amplicon_type %in% c("circular", "possibly_circular")], na.rm = TRUE),
      N = sum(amplicon_type %in% c("circular", "possibly_circular"), na.rm = TRUE)
  ), by = gene_id][!is.na(gene_id) & N > 0]

pdf("plots/TCGA_circos.pdf", width = 10, height = 10)
gcap.plotCircos(TCGA, genome_build = "hg38", 
                highlight_genes = genes_summary[ cn > 40 & N > 1][, .(gene_id, label = N)])
dev.off()

tcga_gene_summary = data.table::copy(genes_summary)
# TCGA Gene cluster ORA analysis -----------------------------------------------

TCGA = readRDS("data/TCGA_SNP.rds")

library(clusterProfiler)
library(gcap)

# ORA_hallmark = gcap.enrich(TCGA, target = "circular")
# ORA_reactome = gcap.enrich(TCGA, target = "circular", category = "C2", subcategory = "CP:REACTOME")
# ORA_KEGG = gcap.enrich(TCGA, target = "circular", category = "C2", subcategory = "CP:KEGG")
# ORA_BP = gcap.enrich(TCGA, target = "circular", category = "C5", subcategory = "GO:BP")
# ORA_MF = gcap.enrich(TCGA, target = "circular", category = "C5", subcategory = "GO:MF")
# ORA_CC = gcap.enrich(TCGA, target = "circular", category = "C5", subcategory = "GO:CC")
# ORA_TFT = gcap.enrich(TCGA, target = "circular", category = "C3", subcategory = "TFT:GTRD")
# ORA_immune = gcap.enrich(TCGA, target = "circular",  category = "C7", subcategory = "IMMUNESIGDB")
# 
# dotplot(ORA_hallmark, showCategory = 20, label_format = 50)
# dotplot(ORA_reactome, showCategory = 20, label_format = 50)
# dotplot(ORA_KEGG, showCategory = 20, label_format = 50)
# dotplot(ORA_BP, showCategory = 20, label_format = 50)
# dotplot(ORA_MF, showCategory = 20, label_format = 50)
# dotplot(ORA_CC, showCategory = 20, label_format = 50)
# dotplot(ORA_TFT, showCategory = 20, label_format = 50)
# dotplot(ORA_immune, showCategory = 20, label_format = 50)
# 
# p = dotplot(ORA_BP, showCategory = 20, label_format = 50)
# ggsave("plots/TCGA_GO_BP.pdf", p, width = 10, height = 9)
# 
# p = dotplot(ORA_reactome, showCategory = 20, label_format = 50)
# ggsave("plots/TCGA_reactome.pdf", p, width = 10, height = 9)
TCGA$sample_summary$case = substr(TCGA$sample_summary$sample, 1, 15)

TCGA$sample_summary = dplyr::left_join(
  TCGA$sample_summary,
  data.table::as.data.table(UCSCXenaShiny::tcga_clinical[, c("sample", "type")]) %>% unique(),
  by = c("case" = "sample")
)

ORA_BP_type_TCGA = gcap.enrich(TCGA, class_by = "type", target = "circular", category = "C5", subcategory = "GO:BP")

ORA_BP_type_TCGA@compareClusterResult = ORA_BP_type_TCGA@compareClusterResult %>% 
  dplyr::mutate(
    Description = stringr::str_remove(Description, "GOBP_") %>% stringr::str_to_title() %>% stringr::str_replace_all("_", " ")
  )

to_show2 = ORA_BP_type_TCGA@compareClusterResult %>%
  dplyr::filter(pvalue < 0.05, p.adjust < 0.2) %>%
  dplyr::count(Description) %>%
  dplyr::filter(n > 2) %>%
  dplyr::pull(Description)

library(ggplot2)
undebug(enrichplot:::dotplot.compareClusterResult)
p = enrichplot:::dotplot.compareClusterResult(ORA_BP_type_TCGA, showCategory = to_show2, label_format = 80) +
  ggpubr::rotate_x_text(45, size = 8) +
  theme(axis.text.y = element_text(size = 8))
p

ggsave("plots/TCGA_BP_by_cancer_type.pdf", p, width = 15, height = 9)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
PCAWG = readRDS("data/PCAWG.rds")
# PCAWG Gene analysis -----------------------------------------------------------

PCAWG$convertGeneID(genome_build = "hg19")
gene_info = PCAWG$gene_summary[!is.na(gene_id)][order(circular, Total, decreasing = TRUE)]
sum(gene_info$circular > 1)

#gcap.plotProfile(PCAWG, top_n = 100, show_column_names = FALSE)
#gcap.plotDistribution(PCAWG, x = gene_info[circular > 2]$gene_id)

gcap.plotCircos(PCAWG, genome_build = "hg19")

genes_summary = PCAWG$data[
  , .(cn = mean(total_cn[amplicon_type %in% c("circular", "possibly_circular")], na.rm = TRUE),
      prob = mean(prob[amplicon_type %in% c("circular", "possibly_circular")], na.rm = TRUE),
      N = sum(amplicon_type %in% c("circular", "possibly_circular"), na.rm = TRUE)
  ), by = gene_id][!is.na(gene_id) & N > 0]

pdf("plots/PCAWG_circos.pdf", width = 10, height = 10)
gcap.plotCircos(PCAWG, genome_build = "hg19", 
                highlight_genes = genes_summary[ cn > 40 & N > 1][, .(gene_id, label = N)])
dev.off()

PCAWG_gene_summary = data.table::copy(genes_summary)
# PCAWG Gene cluster ORA analysis -----------------------------------------------

PCAWG = readRDS("data/PCAWG.rds")

devtools::load_all("~/proj/gcap/")
library(clusterProfiler)

# ORA_hallmark = gcap.enrich(PCAWG, target = "circular")
# ORA_reactome = gcap.enrich(PCAWG, target = "circular", category = "C2", subcategory = "CP:REACTOME")
# ORA_KEGG = gcap.enrich(PCAWG, target = "circular", category = "C2", subcategory = "CP:KEGG")
# ORA_BP = gcap.enrich(PCAWG, target = "circular", category = "C5", subcategory = "GO:BP")
# ORA_MF = gcap.enrich(PCAWG, target = "circular", category = "C5", subcategory = "GO:MF")
# ORA_CC = gcap.enrich(PCAWG, target = "circular", category = "C5", subcategory = "GO:CC")
# ORA_TFT = gcap.enrich(PCAWG, target = "circular", category = "C3", subcategory = "TFT:GTRD")
# ORA_immune = gcap.enrich(PCAWG, target = "circular",  category = "C7", subcategory = "IMMUNESIGDB")
# 
# dotplot(ORA_hallmark, showCategory = 20, label_format = 50)
# dotplot(ORA_reactome, showCategory = 20, label_format = 50)
# dotplot(ORA_KEGG, showCategory = 20, label_format = 50)
# dotplot(ORA_BP, showCategory = 20, label_format = 50)
# dotplot(ORA_MF, showCategory = 20, label_format = 50)
# dotplot(ORA_CC, showCategory = 20, label_format = 50)
# dotplot(ORA_TFT, showCategory = 20, label_format = 50)
# dotplot(ORA_immune, showCategory = 20, label_format = 50)
# 
# p = dotplot(ORA_BP, showCategory = 20, label_format = 50)
# ggsave("plots/PCAWG_GO_BP.pdf", p, width = 10, height = 9)
# 
# p = dotplot(ORA_reactome, showCategory = 20, label_format = 50)
# ggsave("plots/PCAWG_reactome.pdf", p, width = 10, height = 9)

#undebug(gcap.enrich)
ORA_BP_type = gcap.enrich(PCAWG, class_by = "cancer_type", target = "circular", category = "C5", subcategory = "GO:BP")
dotplot(ORA_BP_type, showCategory = 5, label_format = 50) + ggpubr::rotate_x_text(45) +
  theme(axis.text.y = element_text(size = 6))

ORA_BP_type@compareClusterResult = ORA_BP_type@compareClusterResult %>% 
  dplyr::mutate(
    Description = stringr::str_remove(Description, "GOBP_") %>% stringr::str_to_title() %>% stringr::str_replace_all("_", " ")
  )

to_show = ORA_BP_type@compareClusterResult %>%
  dplyr::filter(pvalue < 0.05, p.adjust < 0.2) %>%
  dplyr::count(Description) %>%
  dplyr::filter(n > 2) %>%
  dplyr::pull(Description)

library(ggplot2)
p = enrichplot:::dotplot.compareClusterResult(ORA_BP_type, showCategory = to_show, label_format = 80) +
  ggpubr::rotate_x_text(45, size = 8) +
  theme(axis.text.y = element_text(size = 8))
p

ggsave("plots/PCAWG_BP_by_cancer_type.pdf", p, width = 13, height = 7)

intersect(to_show, to_show2)
# [1] "Collagen catabolic process"             "Intermediate filament based process"    "Intermediate filament organization"    
# [4] "Molting cycle"                          "Regulation of t cell apoptotic process"
# 胶原蛋白分解代谢的过程
# 中间丝
# 中间纤维组织
# 蜕皮周期
# t细胞凋亡过程的调控 (容易导致T细胞耗竭？)

openxlsx::write.xlsx(
  list(TCGA = tcga_gene_summary[order(N,  cn, prob, decreasing = TRUE)],
       PCAWG = PCAWG_gene_summary[order(N,  cn, prob, decreasing = TRUE)]), 
  file = "data/pancan_predicted_circular_gene_summary.xlsx"
)

