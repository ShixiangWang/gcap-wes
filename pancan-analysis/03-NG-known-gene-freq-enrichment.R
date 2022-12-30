# Compare gene set enriched in ecDNA and linear amplicons

# Enrichment analysis -----------------------------------------------------

library(dplyr)
library(clusterProfiler)
library(org.Hs.eg.db)
library(ggplot2)

freq_dt = readRDS("modeling/data/model_gene_amplicon_freq.rds")[type %in% c("Linear", "Circular")]
freq_dt[freq > 0.3113325 * 10]

#freq_dt[gene_id == "ENSG00000124610"]

geneList_linear <- freq_dt[freq > 0.3113325 * 10 & type == "Linear"]$gene_id
geneList_circular <- freq_dt[freq > 0.3113325 * 10 & type == "Circular"]$gene_id

geneList = list(
  linear = unique(bitr(geneList_linear, "ENSEMBL", "ENTREZID", OrgDb = org.Hs.eg.db)$ENTREZID),
  circular = unique(bitr(geneList_circular, "ENSEMBL", "ENTREZID", OrgDb = org.Hs.eg.db)$ENTREZID)
)

ck_go <- compareCluster(geneCluster = geneList, fun = enrichGO, ont = "BP", OrgDb = org.Hs.eg.db, readable = TRUE)
ck_go = simplify(ck_go)
dotplot(ck_go, showCategory = 20)

p = dotplot(ck_go, showCategory = 20)
ggsave("pancan-analysis/plots/circular_and_linear_gene_cluster_GO_BP_comparison.pdf", p, width = 8, height = 10)

View(ck_go@compareClusterResult)

saveRDS(ck_go@compareClusterResult, file = "pancan-analysis/data/circular_and_linear_gene_cluster_GO_BP_comparison.rds")

circular_and_linear_gene_cluster_GO_BP_comparison = ck_go@compareClusterResult
z = circular_and_linear_gene_cluster_GO_BP_comparison %>% 
  group_by(Cluster) %>% 
  summarise(ls = stringr::str_split(geneID, pattern = "/") %>% 
              unlist()) %>% 
  group_split() %>% 
  purrr::map(~unique(.$ls))
sort(setdiff(z[[2]], z[[1]]))
# Part of 
# "H1-1"      "H1-2"     
# [371] "H1-3"      "H1-4"      "H1-5"      "H1-6"      "H2AC1"     "H2AC11"    "H2AC12"    "H2AC13"    "H2AC15"    "H2AC16"   
# [381] "H2AC17"    "H2AC18"    "H2AC19"    "H2AC20"    "H2AC21"    "H2AC4"     "H2AC6"     "H2AC7"     "H2AC8"     "H2AJ"     
# [391] "H2AW"      "H2BC1"     "H2BC10"    "H2BC11"    "H2BC12"    "H2BC13"    "H2BC14"    "H2BC15"    "H2BC17"    "H2BC18"   
# [401] "H2BC21"    "H2BC3"     "H2BC4"     "H2BC5"     "H2BC6"     "H2BC9"     "H2BE1"     "H2BU1"     "H3-3A"     "H3-4"     
# [411] "H3C1"      "H3C10"     "H3C11"     "H3C12"     "H3C13"     "H3C14"     "H3C15"     "H3C2"      "H3C3"      "H3C4"     
# [421] "H3C6"      "H3C7"      "H3C8"      "H4-16"     "H4C1"      "H4C11"     "H4C12"     "H4C13"     "H4C14"     "H4C15"    
# [431] "H4C3"      "H4C4"      "H4C5"      "H4C7"      "H4C8"      "H4C9"      "HAS2"      "HAX1"      "HCK"       "HCST"     
# [441] "HELB"      "HEY1"      "HFE"       "HGF"       "HILPDA"    "HIPK2"     "HLA-A"     "HLA-B"     "HLA-C"     "HLA-DMB"  
# [451] "HLA-DOA"   "HLA-DOB"   "HLA-DPA1"  "HLA-DPB1"  "HLA-DRB1"  "HLA-E"     "HLA-F"     "HLA-G"     "HLX"

# go_bp_linear <- enrichGO(gene = geneList$linear, OrgDb = org.Hs.eg.db, ont = "BP", readable = TRUE)
# go_bp_circular <- enrichGO(gene = geneList$circular, OrgDb = org.Hs.eg.db, ont = "BP", readable = TRUE)
# 
# dotplot(go_bp_linear, showCategory = 10)
# heatplot(go_bp_linear, showCategory = 10)
# 
# dotplot(go_bp_circular, showCategory = 20)
# heatplot(go_bp_circular, showCategory = 20)

# p = heatplot(go_bp_circular, showCategory = 20)
# keep_list = unique(p$data$Gene)
# 
# keep_dt = overlap_dt2[gene_name %in% keep_list & amplicon_classification == "Circular"]
# 
# unique(keep_dt$sample_barcode)
# keep_dt$amplicon_classification = NULL
# 
# df = keep_dt %>% 
#   dplyr::group_by(gene_name) %>% 
#   dplyr::mutate(n = n()) %>% 
#   dplyr::arrange(desc(n), gene_name) %>% 
#   dplyr::ungroup() %>% 
#   dplyr::left_join(
#     overlap_dt[, .(gene_name, chr, start, end)] %>% unique(),
#     by = "gene_name"
#   )
# 
# saveRDS(df, file = "data/circular_enriched_gene_df.rds")
# print(sort(unique(df$gene_name)))

# anno_gene = df %>%                                                                                                                                                                                                                                                     
#   dplyr::select(gene_name) %>% 
#   dplyr::distinct() %>% 
#   dplyr::left_join(annotables::grch37 %>% dplyr::select(symbol, description), by = c("gene_name" = "symbol")) %>% 
#   unique() %>% 
#   dplyr::distinct(gene_name, .keep_all = TRUE)
# 
# unique(sort(df$gene_name))


# ck_kegg <- compareCluster(geneCluster = geneList, fun = enrichKEGG)
# dotplot(ck_kegg, showCategory = NULL)
# cnetplot(ck_kegg)


# library(ReactomePA)
# ck_reactome <- compareCluster(geneCluster = geneList, fun = enrichPathway, readable = TRUE)
# dotplot(ck_reactome, showCategory = 20)
# cnetplot(ck_reactome)


# Circular amplicon specific gene identification --------------------------

amplicon = readRDS("preprocessing-and-EDA/data/amplicon_hg19.rds")
amplicon = amplicon[amplicon_classification %in% c("Linear", "Circular")]

ref = readRDS("preprocessing-and-EDA/data/hg19_gene_info.rds")

ref[, gene_id := substr(gene_id, 1, 15)]
ref = ref[gene_type == "protein_coding"]
ref$gene_type = NULL
ref$strand = NULL

regions = modules::use("lib/regions.R")
overlap_dt = regions$overlaps(amplicon, ref)

overlap_dt = overlap_dt[intersect_ratio >= 1]
overlap_dt2 = overlap_dt[, .(sample_barcode, gene_id, gene_name, amplicon_classification)]

overlap_dt2

tumor_info = readxl::read_excel("preprocessing-and-EDA/data/ICGC-ecDNA.xlsx", sheet = 3) %>%
  dplyr::filter(tumor_or_normal == "tumor") %>%
  dplyr::select(sample_barcode, lineage)

df = tumor_info %>%
  dplyr::left_join(overlap_dt2, by = "sample_barcode") %>% 
  data.table::as.data.table()

data.table::fwrite(df[lineage %in% c("Colorectal", "Gastric", "Esophageal") & amplicon_classification == "Circular"],
                   file = "pancan-analysis/data/NG_known_GI_circular_genes.csv", sep = ",")

rm(overlap_dt, overlap_dt2, tumor_info)

df_freq = data.table::dcast(
  df[!is.na(gene_name), .N, 
     by = .(gene_name, amplicon_classification)],
  gene_name ~ amplicon_classification, fill = 0L)[order(Circular, Linear, decreasing = TRUE)]
df_freq[, diff := Circular - Linear]

df_freq[Circular > 30 & abs(diff) > 10]

p = ggplot(data = df_freq, aes(x = diff, y = Circular)) +
  geom_point(alpha = 0.5, size = 1.2, col = "black") +
  geom_point(data = subset(df_freq, Circular > 30 & diff > 10), alpha = 0.8, size = 1.2, col = "red") +
  geom_point(data = subset(df_freq, Circular > 30 & diff < -10), alpha = 0.6, size = 1.2, col = "blue") +
  ggrepel::geom_label_repel(
    aes(label = gene_name), size = 2,
    data = df_freq[Circular > 30 & abs(diff) > 10], max.overlaps = 100
  ) + 
  labs(x = "Frequency difference of gene on circular and linear DNA", y = "Frequency of gene on circular DNA") +
  theme(plot.title = element_text(hjust = 0.4)) +
  #geom_hline(yintercept = 30, lty = 4, lwd = 0.6, alpha = 0.8) +
  #geom_vline(xintercept = c(-10, 10), lty = 4, lwd = 0.6, alpha = 0.8) +
  theme_bw() +
  theme(
    panel.border = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(colour = "black")
  ) +
  theme_bw(base_size = 15)
p

ggsave("pancan-analysis/plots/circular_and_linear_gene_frequency_comparison.pdf", p, width = 8, height = 6)

df_freq2 = df_freq[Circular > 30 & abs(diff) > 10]
df_freq2

df2 = df[gene_name %in% df_freq2$gene_name & amplicon_classification == "Circular"]
df2 = df2[, .N, by = .(gene_id, gene_name, lineage)][order(N, decreasing = TRUE)]
df2[, type_dist := paste0(lineage, " (", N, ")")]
df2 = df2[, .(type_dist = paste(type_dist, collapse = ",")), by = .(gene_id, gene_name)]

df3 = merge(df_freq2, df2, by = "gene_name")[order(Circular, decreasing = TRUE)]

saveRDS(df3, file = "pancan-analysis/data/high_freq_circular_gene_df.rds")

rm(list = ls())


# Cancer distribution -----------------------------------------------------
# Visualize gene distribution across cancer types.

df = readRDS("pancan-analysis/data/high_freq_circular_gene_df.rds")
genes = unique(df$gene_id)

library(dplyr)
df2 = df %>% 
  tidyr::separate_rows(type_dist, sep = ",") %>% 
  tidyr::separate(type_dist, c("type", "N"), sep = " \\(") %>% 
  dplyr::mutate(Frequency = as.integer(stringr::str_remove(N, "\\)")))

type_order = df2 %>% dplyr::group_by(type) %>% dplyr::summarise(N = length(gene_name)) %>% 
  dplyr::arrange(dplyr::desc(N)) %>% dplyr::pull(type)

gene_order = df2 %>% dplyr::group_by(gene_name) %>% dplyr::summarise(N = length(type)) %>% 
  dplyr::arrange(dplyr::desc(N)) %>% dplyr::pull(gene_name)

library(ggplot2)

p <- ggplot(df2 %>% 
              dplyr::mutate(type = factor(type, levels = type_order),
                            gene_name = factor(gene_name, levels = gene_order)),
            aes(x = gene_name, y = type, fill = Frequency)) +
  geom_tile() +
  ggpubr::theme_pubr() +
  ggpubr::rotate_x_text(45) +
  labs(x = NULL, y = NULL) +
  scale_fill_viridis_b()
p

ggsave("pancan-analysis/plots/high_freq_amplicon_gene_cancer_type_distribution.pdf", p, width = 14, height = 7)

# TODO: need to analyze in cytoband level?
