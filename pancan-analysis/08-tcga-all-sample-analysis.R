library(gcap)
library(data.table)
library(IDConverter)
options(IDConverter.datapath = system.file("extdata", package = "IDConverter"))

setwd(file.path(PROJ_DIR, "pancan-analysis"))

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Merge PCAWG and TCGA samples and convert them to TCGA IDs if available
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

TCGA = readRDS("data/TCGA_SNP.rds")
PCAWG = readRDS("data/PCAWG.rds")

info_tcga = TCGA$sample_summary
info_pcawg = PCAWG$sample_summary

info_pcawg$sample = convert_icgc(info_pcawg$sample, to = "submitted_sample_id")
sum(is.na(info_pcawg$sample))
sum(startsWith(info_pcawg$sample, "TCGA"))

info_pcawg_tcga = info_pcawg[startsWith(info_pcawg$sample, "TCGA")]
length(intersect(substr(info_pcawg_tcga$sample, 1, 15), substr(info_tcga$sample, 1, 15)))
uniqLen((substr(info_pcawg_tcga$sample, 1, 15)))
uniqLen((substr(info_tcga$sample, 1, 15)))

# Intersection analysis
isps = intersect(substr(info_pcawg_tcga$sample, 1, 15), substr(info_tcga$sample, 1, 15))
isps_tcga = info_tcga[substr(sample, 1, 15) %in% isps, ]
isps_pcawg = info_pcawg_tcga[substr(sample, 1, 15) %in% isps, ]
isps_dt = merge(isps_tcga[, .(sample = substr(sample, 1, 15), class)],
                isps_pcawg[, .(sample = substr(sample, 1, 15), class)], by = "sample")
sum(isps_dt$class.x == isps_dt$class.y) / nrow(isps_dt)
#0.7548544

# > nrow(isps_dt)
# [1] 412

info_dt = rbind(
  info_tcga[, .(sample = substr(sample, 1, 15), class, source = "tcga")],
  info_pcawg_tcga[, .(sample = substr(sample, 1, 15), class, source = "pcawg")]
)

saveRDS(info_dt, file = "data/TCGA_combined_class.rds")

info_dt = readRDS("data/TCGA_combined_class.rds")
info_dt[, class := set_default_factor(class)]

info_dt2 = info_dt[, .N, by = .(sample, class, source)][order(N, decreasing = TRUE)]
info_dt2

collapse_many = function(dt) {
  l = kit::uniqLen(dt$class)
  if (l == 1) {
    v = as.character(dt$class[1])
  } else {
    idx = which.max(table(3-as.integer(dt$class)))
    v = c("circular", "noncircular", "nofocal")[idx]
  }
  return(v)
}
info_dt3 = info_dt2[, .(class = collapse_many(.SD)), by = .(sample)]
table(info_dt3$class)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Merge all phenotype data
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
library(UCSCXenaShiny)

cibersort = readr::read_tsv("/data3/wsx_data/TCGA_pancan_immune/TCGA.Kallisto.fullIDs.cibersort.relative.tsv")
cibersort = cibersort %>% 
  dplyr::mutate(
    sample = substr(gsub("\\.", "-", SampleID), 1, 15)
  ) %>% 
  dplyr::select(-SampleID, -CancerType, -P.value, -Correlation, -RMSE) %>% 
  dplyr::select(sample, dplyr::everything())

tcga_clinical

pancan_immune = readr::read_csv("/data3/wsx_data/TCGA_pancan_immune/mmc2.csv") %>% 
  dplyr::select(-c("Eosinophils...41", "Neutrophils...48",
                   "Neutrophils...60", "Eosinophils...61",
                   "OS", "OS Time",
                   "PFI", "PFI Time")) # Remove repeated rows and survival data
pancan_immune = pancan_immune[, 1:32]
colnames(pancan_immune)[1] = "patient"
pancan_immune$`TCGA Study` = NULL

tcga_info = dplyr::left_join(
  tcga_clinical %>% 
    dplyr::mutate(
      stage = stringr::str_match(ajcc_pathologic_tumor_stage, "Stage\\s+(.*?)[ABC]?$")[, 2]
    ) %>% 
    dplyr::select(sample, patient, type, age = age_at_initial_pathologic_diagnosis,
                  gender, race, stage) %>% 
    unique(),
  pancan_immune, by = "patient") %>% 
  dplyr::filter(!is.na(sample)) %>% 
  dplyr::select(-patient) %>% 
  dplyr::full_join(
    cibersort %>% dplyr::distinct(sample, .keep_all = TRUE), by = "sample"
  ) %>% 
  dplyr::full_join(
    info_dt3, by = "sample"
  )

sum(duplicated(tcga_info$sample))
saveRDS(tcga_info, file = "data/tcga_tidy_info.rds")

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 先把变量分为多组，然后分析数据用热图展示
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# [1] "sample"                                  "type"                                   
# [3] "age"                                     "gender"                                 
# [5] "race"                                    "stage"                                  
# [7] "Immune Subtype"                          "TCGA Subtype"                           
# [9] "Leukocyte Fraction"                      "Stromal Fraction"                       
# [11] "Intratumor Heterogeneity"                "TIL Regional Fraction"                  
# [13] "Proliferation"                           "Wound Healing"                          
# [15] "Macrophage Regulation"                   "Lymphocyte Infiltration Signature Score"
# [17] "IFN-gamma Response"                      "TGF-beta Response"                      
# [19] "SNV Neoantigens"                         "Indel Neoantigens"                      
# [21] "Silent Mutation Rate"                    "Nonsilent Mutation Rate"                
# [23] "Number of Segments"                      "Fraction Altered"                       
# [25] "Aneuploidy Score"                        "Homologous Recombination Defects"       
# [27] "BCR Evenness"                            "BCR Shannon"                            
# [29] "BCR Richness"                            "TCR Shannon"                            
# [31] "TCR Richness"                            "TCR Evenness"                           
# [33] "CTA Score"                               "Th1 Cells"                              
# [35] "Th2 Cells"                               "Th17 Cells"                             
# [37] "B.cells.naive"                           "B.cells.memory"                         
# [39] "Plasma.cells"                            "T.cells.CD8"                            
# [41] "T.cells.CD4.naive"                       "T.cells.CD4.memory.resting"             
# [43] "T.cells.CD4.memory.activated"            "T.cells.follicular.helper"              
# [45] "T.cells.regulatory..Tregs."              "T.cells.gamma.delta"                    
# [47] "NK.cells.resting"                        "NK.cells.activated"                     
# [49] "Monocytes"                               "Macrophages.M0"                         
# [51] "Macrophages.M1"                          "Macrophages.M2"                         
# [53] "Dendritic.cells.resting"                 "Dendritic.cells.activated"              
# [55] "Mast.cells.resting"                      "Mast.cells.activated"                   
# [57] "Eosinophils"                             "Neutrophils"                            
# [59] "class"                                  

tcga_info = readRDS(file = "data/tcga_tidy_info.rds")

dt = tcga_info[!is.na(tcga_info$class), ]
dt$Th1_Th2_ratio = scales::rescale(dt$`Th1 Cells`, c(1, 10)) / scales::rescale(dt$`Th2 Cells`, c(1, 10))

var_list = list(
  cibersort = colnames(cibersort)[-1],
  signature = colnames(dt)[c(9:18, 33:36, 60)],
  lesion = colnames(dt)[19:26],
  receptor = colnames(dt)[27:32]
)

rvlist = purrr::map(var_list, function(ls) {
  ge = group_enrichment(dt,
                        grp_vars = "class",
                        enrich_vars = ls,
                        co_method = "wilcox.test")
  ge$grp1 = set_default_factor(ge$grp1)
  ge$enrich_var = factor(ge$enrich_var, ls)
  p = show_group_enrichment(ge, cut_p_value = TRUE, cluster_row = FALSE, return_list = TRUE)
  p = p[[1]] + ggpubr::rotate_x_text(30) + labs(x = NULL, y = NULL)
  return(list(data = ge, plot = p))
})

ggsave("plots/tcga_fCNA_cibersort.pdf", rvlist$cibersort$plot, width = 10, height = 4, device = cairo_pdf)
ggsave("plots/tcga_fCNA_signature.pdf", rvlist$signature$plot, width = 10, height = 4, device = cairo_pdf)
ggsave("plots/tcga_fCNA_lesion.pdf", rvlist$lesion$plot, width = 6, height = 4, device = cairo_pdf)
ggsave("plots/tcga_fCNA_receptor.pdf", rvlist$receptor$plot, width = 5, height = 3, device = cairo_pdf)

rvlist2 = purrr::map(var_list, function(ls) {
  ge = group_enrichment(dt,
                        grp_vars = "class",
                        enrich_vars = ls,
                        co_method = "wilcox.test",
                        ref_group = "noncircular")
  # ge$grp1 = set_default_factor(ge$grp1)
  ge$enrich_var = factor(ge$enrich_var, ls)
  p = show_group_enrichment(ge, cut_p_value = TRUE, cluster_row = FALSE, return_list = TRUE)
  p = p[[1]] + ggpubr::rotate_x_text(30) + labs(x = NULL, y = NULL)
  return(list(data = ge, plot = p))
})

rvlist2$cibersort$plot
rvlist2$signature$plot
rvlist2$lesion$plot
rvlist2$receptor$plot

ggsave("plots/tcga_fCNA_cibersort_noncircular_as_ref.pdf", rvlist2$cibersort$plot, 
       width = 10, height = 3.5, device = cairo_pdf)
ggsave("plots/tcga_fCNA_signature_noncircular_as_ref.pdf", rvlist2$signature$plot,
       width = 10, height = 3.5, device = cairo_pdf)
ggsave("plots/tcga_fCNA_lesion_noncircular_as_ref.pdf", rvlist2$lesion$plot,
       width = 6, height = 3.5, device = cairo_pdf)
ggsave("plots/tcga_fCNA_receptor_noncircular_as_ref.pdf", rvlist2$receptor$plot,
       width = 5, height = 2.7, device = cairo_pdf)


saveRDS(rvlist, file = "data/immune_analysis.rds")
saveRDS(rvlist2, file = "data/immune_analysis_noncircular_as_ref.rds")

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Mutational Signatures Analysis
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
SBS = data.table::fread("/data3/wsx_data/MutationalSignatures/TCGA_WES_sigProfiler_SBS_signatures_in_samples.csv")
ID = data.table::fread("/data3/wsx_data/MutationalSignatures/TCGA_WES_sigProfiler_ID_signatures_in_samples.csv")

dt_SBS = info_tcga %>% 
  dplyr::select(sample, class) %>% 
  dplyr::left_join(SBS %>% 
                     dplyr::mutate(sample = substr(`Sample Names`, 1, 15)) %>% 
                     dplyr::select(-c("Cancer Types", "Accuracy", "Sample Names")) %>% 
                     dplyr::distinct(sample, .keep_all = TRUE), by = "sample") %>% 
  dplyr::filter(!is.na(SBS1))

dt_ID = info_tcga %>% 
  dplyr::select(sample, class) %>% 
  dplyr::left_join(ID %>% 
                     dplyr::mutate(sample = substr(`Sample Names`, 1, 15)) %>% 
                     dplyr::select(-c("Cancer Types", "Accuracy", "Sample Names")) %>% 
                     dplyr::distinct(sample, .keep_all = TRUE), by = "sample") %>% 
  dplyr::filter(!is.na(ID1))

dt_SBS2 = data.table::copy(dt_SBS)
dt_SBS2[, MSI := SBS6 + SBS15 + SBS21 + SBS26 + SBS44]
dt_SBS2[, MSI := ifelse(MSI > 0, "MSI", "MSS")]
dt_SBS2[, class := factor(class, c("nofocal", "noncircular", "circular"))]
msi_rv = sigminer::get_group_comparison(dt_SBS2, 
                                        col_group = "class", cols_to_compare = "MSI", type = "ca")
msi_rv_p = show_group_comparison(msi_rv,
                                 ca_p_threshold = 0.001, 
                                 legend_position_ca = "top", legend_title_ca = "class",
                                 xlab = NA,
                                 set_ca_sig_yaxis = T, text_angle_x = 0, text_hjust_x = 0.5)
p = msi_rv_p$ca$MSI +
  ggplot2::scale_fill_manual(values = c(nofocal = "grey", noncircular = "blue", circular = "red"),
                             name = "class") + labs(x = NULL)
p
ggsave(filename = "../manuscript/plots/TCGA_MSI_by_mutsig.pdf", p, width = 4, height = 4)


call_grp_analysis = function(dt, ref_group = NA) {
  ge = group_enrichment(dt,
                        grp_vars = "class",
                        enrich_vars = colnames(dt)[-c(1, 2)],
                        co_method = "wilcox.test",
                        ref_group = ref_group)
  if (is.na(ref_group)) ge$grp1 = set_default_factor(ge$grp1)
  ge$enrich_var = factor(ge$enrich_var, colnames(dt)[-c(1, 2)])
  p = show_group_enrichment(ge, cut_p_value = TRUE, cluster_row = FALSE, return_list = TRUE)
  p = p[[1]] + labs(x = NULL, y = NULL) + ggpubr::rotate_x_text() #+ coord_flip()
  return(list(data = ge, plot = p))
}

rv_SBS = call_grp_analysis(dt_SBS)
rv_SBS$plot

ggsave("plots/tcga_fCNA_SBS.pdf", rv_SBS$plot,
       width = 22, height = 2.5, device = cairo_pdf)

rv_ID = call_grp_analysis(dt_ID)
rv_ID$plot
ggsave("plots/tcga_fCNA_ID.pdf", rv_ID$plot,
       width = 8, height = 2.5, device = cairo_pdf)

saveRDS(list(SBS = dt_SBS, ID = dt_ID), file = "data/fCNA_MutationalSignature.rds")

###### CN signatures

TCGA = readRDS("data/TCGA_SNP.rds")
info_tcga = TCGA$sample_summary

CNS = info_tcga %>% dplyr::select(sample, dplyr::starts_with("CN", ignore.case = FALSE)) %>% 
  dplyr::mutate(sample = substr(sample, 1, 15)) %>% 
  dplyr::group_by(sample) %>% 
  dplyr::summarise_all(mean, na.rm = TRUE)

dt_CNS = dt %>% 
  dplyr::select(sample, class) %>% 
  dplyr::left_join(CNS, by = "sample")

head(dt_CNS)

rv_CNS = call_grp_analysis(dt_CNS)
rv_CNS$plot
ggsave("plots/tcga_fCNA_CNS.pdf", rv_CNS$plot,
       width = 8, height = 2.5, device = cairo_pdf)

saveRDS(rv_CNS, file = "data/fCNA_CNSignature.rds")


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Immune Related Analysis
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# df_cibersort = data.table::as.data.table(tcga_info)[!is.na(class), c(colnames(cibersort), "class"), with = FALSE]
# df_cibersort = melt(df_cibersort, id.vars = c("sample", "class"), variable.name = "cell_type", value.name = "fraction")
# df_cibersort[, class := set_default_factor(class)]

# library(ggstatsplot)
# library(gginnards)
# p = grouped_ggbetweenstats(
#   data = df_cibersort,
#   x = class,
#   y = fraction,
#   type = "np",
#   centrality.type = "np",
#   grouping.var = cell_type,
#   plotgrid.args = list(nrow = 6),
#   palette = "default_jama",
#   package = "ggsci",
#   bf.message = FALSE,
#   ggplot.component = list(labs(y = NULL)),
#   ggtheme = cowplot::theme_cowplot() + theme(plot.subtitle=element_text(size=7)), 
# )

# library(rstatix)
# stat.test <- df_cibersort %>%
#   group_by(cell_type) %>%
#   t_test(fraction ~ class) %>% 
#   add_y_position()
# stat.test
# 
# library(ggpubr)
# p = ggboxplot(df_cibersort, x = "class", y = "fraction", 
#           fill = "#FC4E07", facet.by = "cell_type", scales = "free", width = 0.5) +
#   theme(legend.position = "none") +
#   stat_pvalue_manual(stat.test, label = "p.adj", tip.length = 0.01) +
#   scale_y_continuous(expand = expansion(mult = c(0.05, 0.1)))
# p
# ggsave("plots/tcga_immune_analysis_cibersort.pdf", p, width = 15, height = 14)
# 
# df1 = data.table::as.data.table(tcga_info)[!is.na(class),
#                                            colnames(tcga_info)[ c(1, 9:26, 33:36, 59)], with =FALSE]
# df1 = melt(df1, id.vars = c("sample", "class"), variable.name = "measure", value.name = "score")
# df1[, class := set_default_factor(class)]
# 
# stat.test <- df1 %>%
#   group_by(measure) %>%
#   t_test(score ~ class) %>% 
#   add_y_position(scales = "free_y", step.increase = 0.2)
# stat.test
# 
# my_comparison = list(
#   c("nofocal", "noncircular"),
#   c("nofocal", "circular"),
#   c("noncircular", "circular")
# )
# p = ggboxplot(df1, x = "class", y = "score", 
#               fill = "#FC4E07", facet.by = "measure", scales = "free", width = 0.5) +
#   # stat_pvalue_manual(stat.test, label = "p.adj.signif") +
#   # scale_y_continuous(expand = expansion(mult = c(0.05, 0.1)))
#   #theme(legend.position = "none") +
#   stat_compare_means(comparisons = my_comparison, method = "wilcox.test", label = "p.signif", tip.length = 0) +
#   scale_y_continuous(expand = expansion(mult = c(0.05, 0.15)))
# p
