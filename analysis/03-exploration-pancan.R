PROJ_DIR = "~/gcap-analysis/manuscript/"
setwd(PROJ_DIR)

library(dplyr)
library(gcaputils)
library(ggplot2)

PCAWG = readRDS("data/PCAWG.rds")

data_tcga = readRDS("data/tcga_snp_data.rds")
data_pcawg = PCAWG$sample_summary
data_tcga$class = set_default_factor(data_tcga$class)
data_pcawg$class = set_default_factor(data_pcawg$class)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# CIN and some overall scores
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
source("../lib/plot.R")

data_tcga$cnv = rowSums(data_tcga[, paste0("CN", 1:19)])

p = plot_clean_stat_box(data_tcga %>% mutate(cnv = log2(cnv + 1e-6)), x = class, y = cnv,
                        ylab = "CN segments (log2 based)")
p

ggsave("plots/TCGA_CNV_events_by_class.pdf", p, width = 5, height = 4)

p = plot_clean_stat_box(data_tcga, x = class, y = cna_burden,
                        ylab = "CNA burden")
p

ggsave("plots/TCGA_CNA_burden_by_class.pdf", p, width = 5, height = 4)

p = plot_clean_stat_box(data_tcga, x = class, y = pLOH,
                        ylab = "Genome percentage with LOH")
p

ggsave("plots/TCGA_pLOH_by_class.pdf", p, width = 5, height = 4)

p = plot_clean_stat_box(data_tcga, x = class, y = AScore,
                        ylab = "Aneuploidy score")
p

ggsave("plots/TCGA_Aneuploidy_score_by_class.pdf", p, width = 5, height = 4)

p = plot_clean_stat_box(data_tcga, x = class, y = purity,
                        ylab = "Tumor purity")
p
ggsave("plots/TCGA_purity_by_class.pdf", p, width = 5, height = 4)

p = plot_clean_stat_box(data_tcga, x = class, y = ploidy,
                        ylab = "Tumor ploidy")
p

ggsave("plots/TCGA_ploidy_by_class.pdf", p, width = 5, height = 4)

# PCAWG
data_pcawg$cnv = rowSums(data_pcawg[, paste0("CN", 1:19)])

p = plot_clean_stat_box(data_pcawg %>% mutate(cnv = log2(cnv + 1e-6)), x = class, y = cnv,
                        ylab = "CN segments (log2 based)")
p
ggsave("plots/PCAWG_CNV_events_by_class.pdf", p, width = 5, height = 4)

p = plot_clean_stat_box(data_pcawg, x = class, y = cna_burden,
                        ylab = "CNA burden")
p

ggsave("plots/PCAWG_CNA_burden_by_class.pdf", p, width = 5, height = 4)

p = plot_clean_stat_box(data_pcawg, x = class, y = pLOH,
                        ylab = "Genome percentage with LOH")
p

ggsave("plots/PCAWG_pLOH_by_class.pdf", p, width = 5, height = 4)

p = plot_clean_stat_box(data_pcawg, x = class, y = AScore,
                        ylab = "Aneuploidy score")
p

ggsave("plots/PCAWG_Aneuploidy_score_by_class.pdf", p, width = 5, height = 4)


p = plot_clean_stat_box(data_pcawg, x = class, y = purity,
                        ylab = "Tumor purity")
p

ggsave("plots/PCAWG_purity_by_class.pdf", p, width = 5, height = 4)

p = plot_clean_stat_box(data_pcawg, x = class, y = ploidy,
                        ylab = "Tumor ploidy")
p

ggsave("plots/PCAWG_ploidy_by_class.pdf", p, width = 5, height = 4)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Immune landscape and others
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
library(gcaputils)
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
    dplyr::select(sample, patient, age = age_at_initial_pathologic_diagnosis,
                  gender, race, stage) %>% 
    unique(),
  pancan_immune, by = "patient") %>% 
  dplyr::filter(!is.na(sample)) %>% 
  dplyr::select(-patient) %>% 
  dplyr::full_join(
    cibersort %>% dplyr::distinct(sample, .keep_all = TRUE), by = "sample"
  ) %>% 
  dplyr::right_join(
    data_tcga %>% dplyr::select(-c(OS:cnv), type), by = "sample"
  )

sum(duplicated(tcga_info$sample))
saveRDS(tcga_info, file = "data/tcga_tidy_info.rds")

# 先把变量分为多组，然后分析数据用热图展示

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
# [59] "purity"                                  "ploidy"                                 
# [61] "AScore"                                  "pLOH"                                   
# [63] "cna_burden"                              "CN1"                                    
# [65] "CN2"                                     "CN3"                                    
# [67] "CN4"                                     "CN5"                                    
# [69] "CN6"                                     "CN7"                                    
# [71] "CN8"                                     "CN9"                                    
# [73] "CN10"                                    "CN11"                                   
# [75] "CN12"                                    "CN13"                                   
# [77] "CN14"                                    "CN15"                                   
# [79] "CN16"                                    "CN17"                                   
# [81] "CN18"                                    "CN19"                                   
# [83] "class"                                                 

dt = readRDS(file = "data/tcga_tidy_info.rds")
dt$Th1_Th2_ratio = scales::rescale(dt$`Th1 Cells`, c(1, 10)) / scales::rescale(dt$`Th2 Cells`, c(1, 10))
table(substr(dt$sample, 14, 15))
# 01   02   03   06 
# 9244    3  113  339 
# https://gdc.cancer.gov/resources-tcga-users/tcga-code-tables/sample-type-codes

var_list = list(
  cibersort = colnames(cibersort)[-1],
  signature = colnames(dt)[c(8:17, 32:35, 84)],
  lesion = colnames(dt)[18:25],
  receptor = colnames(dt)[26:31]
)

library(sigminer)
rvlist = purrr::map(var_list, function(ls) {
  ge = group_enrichment(dt,
                        grp_vars = "class",
                        enrich_vars = ls,
                        co_method = "wilcox.test")
  ge$grp1 = set_default_factor(ge$grp1)
  ge$enrich_var = factor(ge$enrich_var, ls)
  ge$fdr = ifelse(ge$fdr == 0, ge$fdr + .Machine$double.xmin, ge$fdr)
  p = show_group_enrichment(ge, cut_p_value = TRUE, cluster_row = FALSE, return_list = TRUE)
  p = p[[1]] + ggpubr::rotate_x_text(30) + labs(x = NULL, y = NULL)
  return(list(data = ge, plot = p))
})

ggsave("plots/tcga_fCNA_cibersort.pdf", rvlist$cibersort$plot, width = 10, height = 3, device = cairo_pdf)
ggsave("plots/tcga_fCNA_signature.pdf", rvlist$signature$plot, width = 10, height = 4, device = cairo_pdf)
ggsave("plots/tcga_fCNA_lesion.pdf", rvlist$lesion$plot, width = 6, height = 3, device = cairo_pdf)
ggsave("plots/tcga_fCNA_receptor.pdf", rvlist$receptor$plot, width = 5, height = 3, device = cairo_pdf)

#undebug(show_group_enrichment)
rvlist2 = purrr::map(var_list, function(ls) {
  ge = group_enrichment(dt,
                        grp_vars = "class",
                        enrich_vars = ls,
                        co_method = "wilcox.test",
                        ref_group = "noncircular")
  # ge$grp1 = set_default_factor(ge$grp1)
  ge$enrich_var = factor(ge$enrich_var, ls)
  ge$fdr = ifelse(ge$fdr == 0, ge$fdr + .Machine$double.xmin, ge$fdr)
  p = show_group_enrichment(ge, cut_p_value = TRUE, cluster_row = FALSE, return_list = TRUE)
  p = p[[1]] + ggpubr::rotate_x_text(30) + labs(x = NULL, y = NULL)
  return(list(data = ge, plot = p))
})

rvlist2$cibersort$plot
rvlist2$signature$plot
rvlist2$lesion$plot
rvlist2$receptor$plot

ggsave("plots/tcga_fCNA_cibersort_noncircular_as_ref.pdf", rvlist2$cibersort$plot, 
       width = 10, height = 3, device = cairo_pdf)
ggsave("plots/tcga_fCNA_signature_noncircular_as_ref.pdf", rvlist2$signature$plot,
       width = 10, height = 3, device = cairo_pdf)
ggsave("plots/tcga_fCNA_lesion_noncircular_as_ref.pdf", rvlist2$lesion$plot,
       width = 6, height = 3, device = cairo_pdf)
ggsave("plots/tcga_fCNA_receptor_noncircular_as_ref.pdf", rvlist2$receptor$plot,
       width = 5, height = 2.7, device = cairo_pdf)


saveRDS(rvlist, file = "data/immune_analysis.rds")
saveRDS(rvlist2, file = "data/immune_analysis_noncircular_as_ref.rds")


dir.create("plots/TCGA_cibersort")
for (i in colnames(dt)[36:57]) {
  dd = dt %>% dplyr::select(c("class", i))
  th = mean(dd[[2]], na.rm = TRUE) + 2*sd(dd[[2]], na.rm = TRUE)
  #th = quantile(dd[[2]], na.rm = TRUE, 0.75) + 1.5*IQR(dd[[2]], na.rm = TRUE)
  dd[[2]] = ifelse(dd[[2]] > th, th, dd[[2]])
  p = plot_clean_stat_box(dd, x = class, y = !!i,
                          ylab = i)
  p
  ggsave(sprintf("plots/TCGA_cibersort/%s_by_class.pdf", i), p, width = 5, height = 4)
}

## TCGA_cells_by_Kassandra
## 这个得到的分类完全没有cibersort精细
KA_dt = lapply(list.files("data/TCGA_cells_by_Kassandra/", full.names = TRUE),
               function(x) {
                 d = data.table::fread(x, data.table = FALSE)
                 rownames(d) = d$V1
                 d$V1 = NULL
                 d = t(d)  
                 data.table::as.data.table(d, keep.rownames = "sample")
               })
KA_dt = rbindlist(KA_dt)
KA_dt = KA_dt[, colnames(KA_dt)[!endsWith(colnames(KA_dt), "_std")], with = FALSE]
saveRDS(KA_dt, file = "data/Kassandra_tidy.rds")

KA_dt =readRDS("data/Kassandra_tidy.rds")
dt2 = merge(data.table(dt)[, list(sample, class, type)], KA_dt, by = "sample")

dir.create("plots/TCGA_Kassandra")
for (i in colnames(dt2)[-c(1:3)]) {
  dd = dt2 %>% dplyr::select(c("class", i))
  th = mean(dd[[2]], na.rm = TRUE) + 2*sd(dd[[2]], na.rm = TRUE)
  dd[[2]] = ifelse(dd[[2]] > th, th, dd[[2]])
  p = plot_clean_stat_box(dd, x = class, y = !!i,
                          ylab = i)
  p
  ggsave(sprintf("plots/TCGA_Kassandra/%s_by_class.pdf", i), p, width = 5, height = 4)
}

dt2_long = melt(dt2,
                id.vars = c("sample", "type", "class"), 
                variable.name = "TME", value.name = "percentages")

library(sigminer)
rvlist = purrr::map(list(TME = colnames(dt2)[-c(1:3)]), function(ls) {
  ge = group_enrichment(dt2,
                        grp_vars = "class",
                        enrich_vars = ls,
                        co_method = "wilcox.test")
  ge$grp1 = set_default_factor(ge$grp1)
  ge$enrich_var = factor(ge$enrich_var, ls)
  ge$fdr = ifelse(ge$fdr == 0, ge$fdr + .Machine$double.xmin, ge$fdr)
  p = show_group_enrichment(ge, cut_p_value = TRUE, cluster_row = FALSE, return_list = TRUE)
  p = p[[1]] + ggpubr::rotate_x_text(30) + labs(x = NULL, y = NULL)
  return(list(data = ge, plot = p))
})
rvlist$TME$plot
ggsave("plots/tcga_fCNA_Kassandra_TME.pdf", rvlist$TME$plot, width = 10, height = 3, device = cairo_pdf)

corrplot::corrplot(cor(dt2[, -c(1:3)]))

# fCNA_rv = data.table()
# for (i in colnames(dt2)[-c(1:3)]) {
#   data = data.table::as.data.table(dt2)[class != "nofocal"][!is.na(class) & !is.na(type)]
#   data = data[type %in% data[, list(N = min(sum(class == "circular"), sum(class != "circular"))), by = list(type)][N >= 3]$type]
#   rv = group_enrichment2(data, 
#                          subset_var = "type", 
#                          grp_vars = "class", 
#                          enrich_vars = i,
#                          co_method = "wilcox.test",
#                          ref_group = "noncircular")
#   rv$grp_var = i
#   fCNA_rv = rbind(fCNA_rv, rv)
# }
# 
# fCNA_rv2 = fCNA_rv[grp1 == "circular"]
# fCNA_rv2[, `:=`(target = grp1, grp1 = grp_var,
#                 grp_var = "circular",
#                 fdr = p.adjust(p_value, "fdr"))]
# p = show_group_enrichment(fCNA_rv2,
#                           cluster_row = FALSE, return_list = TRUE, cut_p_value = TRUE, use_fdr = FALSE)
# p = p$circular + labs(x = NULL, y = NULL)
# p


# library(ggplot2)
# library(jjPlot)
# library(ggnewscale)
# dt2_long_grp = dt2_long[, list(
#   frac_circ_noncirc = mean(percentages[class == "circular"], na.rm = TRUE) / mean(percentages[class == "noncircular"], na.rm = TRUE),
#   frac_noncirc_nofocal = mean(percentages[class == "noncircular"], na.rm = TRUE) / mean(percentages[class == "nofocal"], na.rm = TRUE)
# ), by = list(type, TME)]
# 
# dt2_long_grp[, `:=`(
#   frac_circ_noncirc = ifelse(frac_circ_noncirc > 2, 2, frac_circ_noncirc),
#   frac_noncirc_nofocal = ifelse(frac_noncirc_nofocal > 2, 2, frac_noncirc_nofocal)
# )]
# ggplot(dt2_long_grp,
#        aes(x = type, y = TME)) +
#   geom_jjtriangle(aes(fill = frac_noncirc_nofocal),type = 'ul') +
#   scale_fill_gradient2(low = 'blue', mid = "white", high = 'red',
#                        midpoint = 1) +
#   # new legend
#   new_scale_fill() +
#   geom_jjtriangle(aes(fill = frac_circ_noncirc),type = 'br') +
#   scale_fill_gradient2(low = 'blue', mid = "white", high = '#FF3399', 
#                        midpoint = 1) +
#   coord_fixed() +
#   ggpubr::rotate_x_text()



source("../lib/plot.R")

# plot_clean_stat_box(
#   dt2,
#   x = class, y = Immune_general, ylab = "Leukocyte Fraction"
# )

plot_clean_stat_box(
  dt2 %>% dplyr::mutate(Tregs = ifelse(Tregs > 0.02, 0.02, Tregs)),
  x = class, y = Tregs, ylab = "Tregs"
) -> p
ggsave("plots/TCGA_fCNA_class_Tregs.pdf", p, width = 5, height = 4)

p = plot_clean_stat_box(
  dt,
  x = class, y = `Leukocyte Fraction`, ylab = "Leukocyte Fraction"
)
ggsave("plots/TCGA_fCNA_class_LF.pdf", p, width = 5, height = 4)

p = plot_clean_stat_box(
  dt,
  x = class, y = Proliferation, ylab = "Proliferation"
)
ggsave("plots/TCGA_fCNA_class_Proliferation.pdf", p, width = 5, height = 4)

p = plot_clean_stat_box(
  dt,
  x = class, y = `IFN-gamma Response`, ylab = "IFN-gamma Response"
)
ggsave("plots/TCGA_fCNA_class_IFN_gamma.pdf", p, width = 5, height = 4)

p = plot_clean_stat_box(
  dt,
  x = class, y = `TGF-beta Response`, ylab = "TGF-beta Response"
)
ggsave("plots/TCGA_fCNA_class_TGF_beta.pdf", p, width = 5, height = 4)

p = plot_clean_stat_box(
  dt,
  x = class, y = `Lymphocyte Infiltration Signature Score`, ylab = "Lymphocyte Infiltration Signature Score"
)
ggsave("plots/TCGA_fCNA_class_LF_signature.pdf", p, width = 5, height = 4)

p = plot_clean_stat_box(
  dt,
  x = class, y = `Macrophage Regulation`, ylab = "Macrophage Regulation"
)
ggsave("plots/TCGA_fCNA_class_Macrophage_Regulation.pdf", p, width = 5, height = 4)


p = plot_clean_stat_box(
  dt,
  x = class, y = `Intratumor Heterogeneity`, ylab = "Intratumor Heterogeneity"
)
ggsave("plots/TCGA_fCNA_class_Intratumor_Heterogeneity.pdf", p, width = 5, height = 4)

p = plot_clean_stat_box(
  dt,
  x = class, y = `Th1 Cells`, ylab = "Th1 Cells"
)
ggsave("plots/TCGA_fCNA_class_Th1.pdf", p, width = 5, height = 4)

p = plot_clean_stat_box(
  dt,
  x = class, y = `Th2 Cells`, ylab = "Th2 Cells"
)
p
ggsave("plots/TCGA_fCNA_class_Th2.pdf", p, width = 5, height = 4)

p = plot_clean_stat_box(
  dt,
  x = class, y = `Th17 Cells`, ylab = "Th17 Cells"
)
p
ggsave("plots/TCGA_fCNA_class_Th17.pdf", p, width = 5, height = 4)

p = plot_clean_stat_box(
  dt,
  x = class, y = Th1_Th2_ratio, ylab = "Ratio of Th1 and Th2"
)
p
ggsave("plots/TCGA_fCNA_class_Th1_2_ratio.pdf", p, width = 5, height = 4)

p = plot_clean_stat_box(
  dt %>% dplyr::mutate(`SNV Neoantigens` = log10(`SNV Neoantigens` + 1)),
  x = class, y = `SNV Neoantigens`, ylab = "log10(SNV Neoantigens + 1)"
)
p
ggsave("plots/TCGA_fCNA_class_SNV_neoantigens.pdf", p, width = 5, height = 4)

p = plot_clean_stat_box(
  dt %>% dplyr::mutate(`Indel Neoantigens` = log10(`Indel Neoantigens` + 1)),
  x = class, y = `Indel Neoantigens`, ylab = "log10(Indel Neoantigens + 1)"
)
p
ggsave("plots/TCGA_fCNA_class_Indel_neoantigens.pdf", p, width = 5, height = 4)


# plot_clean_stat_box(
#   dt2 %>% dplyr::mutate(class = ifelse(as.character(dt2$class) %in% "circular", "ecDNA+", "ecDNA-")) %>%
#     dplyr::select(class, Lymphocytes),
#   x = class, y = Lymphocytes
# )


# Go to the support file apm.R
# # TCGA pancan immune subtypes
# grp_ca = get_group_comparison(dt, col_group = "Immune Subtype",
#                               cols_to_compare = "class")
# plist = show_group_comparison(grp_ca, ca_p_threshold = 0.001, 
#                               legend_position_ca = "top", legend_title_ca = "class",
#                               xlab = NA,
#                               set_ca_sig_yaxis = T, text_angle_x = 0, text_hjust_x = 0.5)
# p = plist$ca$class + ggsci::scale_fill_lancet()


# TIDE
fl = list.files("/data3/wsx_data/TCGA_pancan_immune/TIDE_Results/Tumor_Dysf_Excl_scores/", pattern = "TCGA",
                full.names = TRUE)
fl = fl[grepl("OS_base", fl) & grepl("RNASeq", fl)]
tide = rbindlist(lapply(fl, data.table::fread))
tide = unique(tide)
sum(duplicated(tide$V1))

colnames(tide)[1] = "patient"

dt2 = dt %>% dplyr::select(sample, type, class)#, T.cells.CD8, NK.cells.activated)
dt2
dt2 = dplyr::left_join(dt2 %>% dplyr::mutate(patient = substr(sample, 1, 12)),
                 tide, by = "patient")

# apm_dt = readRDS("data/TCGA_MHC_GSVA_pathway_with_fCNA_class.rds")
# dt3 = merge(data.table::as.data.table(dt2 %>% dplyr::select(-patient)), 
#             apm_dt[, list(sample, GO_MHC_I, GO_MHC_II, REACTOME_MHC_I, REACTOME_MHC_II)], by="sample")
# dt3

# library(ggpubr)
# ggscatter(dt3, x = "Dysfunction", y = "T.cells.CD8",
#           add = "reg.line",                         # Add regression line
#           conf.int = TRUE,                          # Add confidence interval
#           color = "class", palette = c("grey", "blue", "red"),           # Color by groups "cyl"
#           shape = "class"                             # Change point shape by groups "cyl"
# )+
#   stat_cor(aes(color = class))           # Add correlation coefficient


sum(!is.na(dt2$Dysfunction))
source("../lib/plot.R")
p = plot_clean_stat_box(dt2, class, Exclusion, ylab = "T cell exclusion")
p
ggsave("plots/TCGA_T_cell_exclusion_by_class.pdf", p, width = 5, height = 4)

p = plot_clean_stat_box(dt2, class, Dysfunction, ylab = "T cell dysfunction")
p
ggsave("plots/TCGA_T_cell_dysfunction_by_class.pdf", p, width = 5, height = 4)

# https://zhuanlan.zhihu.com/p/496992477
p = plot_clean_stat_box(dt2, class, MDSC, ylab = "MDSC")
p
ggsave("plots/TCGA_MDSC_by_class.pdf", p, width = 5, height = 4)

# https://zhuanlan.zhihu.com/p/453842393
p = plot_clean_stat_box(dt2, class, CAF, ylab = "CAF")
p
ggsave("plots/TCGA_CAF_by_class.pdf", p, width = 5, height = 4)

p = plot_clean_stat_box(dt2, class, M2, ylab = "M2")
p
ggsave("plots/TCGA_M2_by_class.pdf", p, width = 5, height = 4)

# 分瘤种
# https://www.nature.com/articles/s41591-018-0136-1
library(sigminer)

fCNA_rv = data.table()
for (i in c("Exclusion", "Dysfunction", "MDSC", "CAF", "M2")) {
  data = data.table::as.data.table(dt2)[class != "nofocal"][!is.na(class) & !is.na(type)]
  data = data[type %in% data[, list(N = min(sum(class == "circular"), sum(class != "circular"))), by = list(type)][N >= 3]$type]
  rv = group_enrichment2(data, 
                         subset_var = "type", 
                         grp_vars = "class", 
                         enrich_vars = i,
                         co_method = "wilcox.test")
  rv$grp_var = i
  fCNA_rv = rbind(fCNA_rv, rv)
}

fCNA_rv2 = fCNA_rv[grp1 == "circular"]
fCNA_rv2[, `:=`(target = grp1, grp1 = grp_var,
                grp_var = "circular",
                fdr = p.adjust(p_value, "fdr"))]
fCNA_rv2[, grp1 := factor(grp1, rev(c("Dysfunction", "Exclusion", "MDSC", "CAF", "M2")))]
p = show_group_enrichment(fCNA_rv2,
                          cluster_row = FALSE, return_list = TRUE, use_fdr = FALSE)
p = p$circular + labs(x = NULL, y = NULL)
p
ggsave("plots/TCGA_fCNA_class_TIDE_circular_vs_noncircular.pdf",
       p, width = 15, height = 4)
#
fCNA_rv = data.table()
for (i in c("Exclusion", "Dysfunction", "MDSC", "CAF", "M2")) {
  data = data.table::as.data.table(dt2)[class != "noncircular"][!is.na(class) & !is.na(type)]
  data = data[type %in% data[, list(N = min(sum(class == "circular"), sum(class != "circular"))), by = list(type)][N >= 3]$type]
  rv = group_enrichment2(data, 
                         subset_var = "type", 
                         grp_vars = "class", 
                         enrich_vars = i,
                         co_method = "wilcox.test")
  rv$grp_var = i
  fCNA_rv = rbind(fCNA_rv, rv)
}

fCNA_rv2 = fCNA_rv[grp1 == "circular"]
fCNA_rv2[, `:=`(target = grp1, grp1 = grp_var,
                grp_var = "circular",
                fdr = p.adjust(p_value, "fdr"))]
fCNA_rv2[, grp1 := factor(grp1, rev(c("Dysfunction", "Exclusion", "MDSC", "CAF", "M2")))]
p = show_group_enrichment(fCNA_rv2,
                          cluster_row = FALSE, return_list = TRUE, use_fdr = FALSE)
p = p$circular + labs(x = NULL, y = NULL)
p
ggsave("plots/TCGA_fCNA_class_TIDE_circular_vs_nofocal.pdf",
       p, width = 15, height = 4)

#
fCNA_rv = data.table()
for (i in c("Exclusion", "Dysfunction", "MDSC", "CAF", "M2")) {
  data = data.table::as.data.table(dt2)[class != "circular"][!is.na(class) & !is.na(type)]
  data = data[type %in% data[, list(N = min(sum(class == "noncircular"), sum(class != "noncircular"))), by = list(type)][N >= 3]$type]
  rv = group_enrichment2(data, 
                         subset_var = "type", 
                         grp_vars = "class", 
                         enrich_vars = i,
                         co_method = "wilcox.test")
  rv$grp_var = i
  fCNA_rv = rbind(fCNA_rv, rv)
}

fCNA_rv2 = fCNA_rv[grp1 == "noncircular"]
fCNA_rv2[, `:=`(target = grp1, grp1 = grp_var,
                grp_var = "noncircular",
                fdr = p.adjust(p_value, "fdr"))]
fCNA_rv2[, grp1 := factor(grp1, rev(c("Dysfunction", "Exclusion", "MDSC", "CAF", "M2")))]
p = show_group_enrichment(fCNA_rv2,
                          cluster_row = FALSE, return_list = TRUE, use_fdr = FALSE)
p = p$noncircular + labs(x = NULL, y = NULL)
p
ggsave("plots/TCGA_fCNA_class_TIDE_noncircular_vs_nofocal.pdf",
       p, width = 15, height = 4)


cor(dt2$Exclusion, dt2$MDSC, use = "pair")
cor(dt2$Exclusion, dt2$CAF, use = "pair")
cor(dt2$Exclusion, dt2$M2, use = "pair")

# Gene level vis ----------------------------------------------------------

TCGA_SNP = readRDS("data/TCGA_SNP.rds")

## Genome heatmap ----------

source("../lib/hl.R")


# CD274
# data[gene_id == "ENSG00000120217" & gene_class == "circular"]


# dt_band_freq[freq >= 5]
# dt_band_freq[, `:=`(
#   chunk = as.integer(sub("(chr[^:]+:[pq])", "", band2)),
#   chr = sub("(chr[^:]+):.*", "\\1", band2)
# )]

# zz = dt_band_freq %>% 
#   dplyr::group_nest(chr) %>% 
#   dplyr::mutate(data = purrr::map(data, function(df) {
#     
#   }))

#debug(gcap.plotGenomeHeatmap)
devtools::load_all("~/proj/gcaputils/")
# pdf("plots/TCGA_genome_heatmap.pdf", width = 9, height = 8)
# gcap.plotGenomeHeatmap(TCGA_SNP, group = "type", bycol = TRUE, sort = TRUE,
#                        highlights = get_highlights(TCGA_SNP), fontsize = 6)
# dev.off()
# 
pdf("plots/TCGA_genome_heatmap_oncogene.pdf", width = 9, height = 8)
#undebug(gcap.plotGenomeHeatmap)
gcap.plotGenomeHeatmap(TCGA_SNP, group = "type", bycol = TRUE, sort = TRUE,
                       highlights = get_highlights(TCGA_SNP,
                                                   gene_freq = 15,
                                                   grp_freq = 2,
                                                   only_oncogenes = TRUE),
                       fontsize = 6)
dev.off()

pdf("plots/TCGA_genome_heatmap_by_row_oncogene.pdf", width = 9, height = 8)
gcap.plotGenomeHeatmap(TCGA_SNP, group = "type", bycol = FALSE, sort = TRUE,
                       highlights = get_highlights(TCGA_SNP,
                                                   gene_freq = 15,
                                                   grp_freq = 2,
                                                   only_oncogenes = TRUE), 
                       fontsize = 6)
dev.off()

PCAWG = readRDS("data/PCAWG.rds")

# PCAWG为TCGA的1/3左右，所以gene freq也降为1/3
pdf("plots/PCAWG_genome_heatmap_by_row_oncogene.pdf", width = 9.5, height = 8)
gcap.plotGenomeHeatmap(PCAWG, group = "cancer_type", bycol = FALSE, sort = TRUE,
                       highlights = get_highlights(PCAWG, group = "cancer_type",
                                                   gene_freq = 5,
                                                   grp_freq = 2, only_oncogenes = TRUE), 
                       fontsize = 6)
dev.off()

# pdf("plots/PCAWG_genome_heatmap.pdf", width = 9, height = 8)
# gcap.plotGenomeHeatmap(PCAWG, group = "cancer_type", bycol = TRUE, sort = TRUE,
#                        highlights = get_highlights(PCAWG, group = "cancer_type", 
#                                                    gene_freq = 10, grp_freq = 2, min_mean_cn = 10),
#                        fontsize = 6)
# dev.off()
# 
pdf("plots/PCAWG_genome_heatmap_oncogene.pdf", width = 9, height = 8)
gcap.plotGenomeHeatmap(PCAWG, group = "cancer_type", bycol = TRUE, sort = TRUE,
                       highlights = get_highlights(PCAWG, group = "cancer_type",
                                                   gene_freq = 5,
                                                   grp_freq = 2,
                                                   only_oncogenes = TRUE),
                       fontsize = 6)
dev.off()

## Oncogene distribution
devtools::load_all("~/proj/gcaputils/")

TCGA_SNP = readRDS("data/TCGA_SNP.rds")
TCGA_SNP$convertGeneID()
# p = gcaputils::gcap.plotDistribution(TCGA_SNP, by_gene = TRUE, genelist = hl_tcga,
#                                  palette = c("#0066CC", "#CC0033"), x_size = 5, bar_width = 1)
# ggsave("plots/TCGA_oncogene_fCNA_dist.pdf", p, width = 18, height = 4)

# Too many, just output a summary table
tcga_oncogene = TCGA_SNP$getGeneSummary()[gene_id %in% gcap::oncogenes$OncogeneName]
tcga_oncogene$Total = NULL
tcga_oncogene$noncircular_T = nrow(TCGA_SNP$sample_summary[class == "noncircular"])
tcga_oncogene$circular_T = nrow(TCGA_SNP$sample_summary[class == "circular"])
tcga_ref = TCGA_SNP$getGeneSummary()[, list(f = mean(circular / noncircular, na.rm = TRUE))]$f
tcga_oncogene = tcga_oncogene[, {
  tbl = matrix(c(circular, noncircular, circular_T - circular, noncircular_T - noncircular), nrow = 2)
  r = fisher.test(tbl, or = tcga_ref)
  #print(r)
  list(noncircular, circular, p = r$p, or = r$estimate %>% as.numeric())
}, by = list(gene_id)]

tcga_oncogene_sum = tcga_oncogene[
  , list(gene_id, p, or, noncirc_freq = noncircular / nrow(TCGA_SNP$sample_summary),
         circ_freq = circular / nrow(TCGA_SNP$sample_summary))
][order(-or)]

tcga_oncogene_sum[, p.adj := p.adjust(p, method = "fdr")]


PCAWG = readRDS("data/PCAWG.rds")
PCAWG$convertGeneID()
# p = gcaputils::gcap.plotDistribution(PCAWG, by_gene = TRUE, genelist = hl_tcga,
#                                      palette = c("#0066CC", "#CC0033"), x_size = 5, bar_width = 1)
# ggsave("plots/PCAWG_oncogene_fCNA_dist.pdf", p, width = 18, height = 4)
# 
# save(hl_tcga, hl_pcawg, file = "data/TCGA_PCAWG_amplifed_oncogene_list.RData")

pcawg_oncogene = PCAWG$getGeneSummary()[gene_id %in% gcap::oncogenes$OncogeneName]
pcawg_oncogene$Total = NULL
pcawg_oncogene$noncircular_T = nrow(PCAWG$sample_summary[class == "noncircular"])
pcawg_oncogene$circular_T = nrow(PCAWG$sample_summary[class == "circular"])
pcawg_oncogene = pcawg_oncogene[, {
  tbl = matrix(c(circular, noncircular, circular_T - circular, noncircular_T - noncircular), nrow = 2)
  r = fisher.test(tbl, or = tcga_ref)
  #print(r)
  list(noncircular, circular, p = r$p, or = r$estimate %>% as.numeric())
}, by = list(gene_id)]

pcawg_oncogene_sum = pcawg_oncogene[
  , list(gene_id, p, or, noncirc_freq = noncircular / nrow(PCAWG$sample_summary),
         circ_freq = circular / nrow(PCAWG$sample_summary))
][order(-or)]

pcawg_oncogene_sum[, p.adj := p.adjust(p, method = "fdr")]
pcawg_oncogene_sum

openxlsx::write.xlsx(
  list(TCGA = tcga_oncogene_sum,
       PCAWG = pcawg_oncogene_sum),
  file = "data/circular_noncircular_tendency_analysis.xlsx"
)

overlap_up = intersect(
  tcga_oncogene_sum[or > tcga_ref & p.adj < 0.05]$gene_id,
  pcawg_oncogene_sum[or > tcga_ref & p.adj < 0.05]$gene_id
)

overlap_down = intersect(
  tcga_oncogene_sum[or < tcga_ref & p.adj < 0.05]$gene_id,
  pcawg_oncogene_sum[or < tcga_ref & p.adj < 0.05]$gene_id
)

saveRDS(c(overlap_down, overlap_up), file = "circular_noncircular_tendency_genes.rds")

# 基本都倾向于noncircular，意义不大

# hl_tcga = unique(get_highlights(readRDS("data/TCGA_SNP.rds"),
#                                 gene_freq = 15,
#                                 grp_freq = 2,
#                                 only_oncogenes = TRUE, return_dt = TRUE)$gene_name)
# hl_pcawg = unique(get_highlights(readRDS("data/PCAWG.rds"),
#                                 gene_freq = 5,
#                                 grp_freq = 2,
#                                 only_oncogenes = TRUE, return_dt = TRUE)$gene_name)


tcga_oncogene_sum = readxl::read_excel("data/circular_noncircular_tendency_analysis.xlsx", sheet = 1)
pcawg_oncogene_sum = readxl::read_excel("data/circular_noncircular_tendency_analysis.xlsx", sheet = 2)

hl = c("ERBB2", "EGFR", "CDK4", "MYC", "FGFR2", "MDM2", "CCND1") 
hl = c("MYC", "MYCL", "MYCN", "EGFR", "ERBB2", "FGFR1", "MET", "MDM2", "MDM4", "CCND1",
       "CCNE1", "SOX2", "E2F3", "CDK6", "KRAS")
# MYC, MYCL, MYCN, EGFR, ERBB2, FGFR1, MET, MDM2, MDM4, CCND1, CCNE1, SOX2, E2F3, CDK6, KRAS


p1 = ggplot(tcga_oncogene_sum,
       aes(x = or, y = -log10(p.adj))) +
  geom_point(size = 0.5, alpha = 0.1) +
  geom_point(data = subset(tcga_oncogene_sum, gene_id %in% hl), color = "red") +
  ggrepel::geom_text_repel(aes(label = gene_id),
                           data = subset(tcga_oncogene_sum, gene_id %in% hl),
                           max.overlaps = 100) +
  geom_vline(xintercept = tcga_ref, linetype = 2) +
  cowplot::theme_cowplot() +
  labs(x = "Odds ratio")
p1

p2 = ggplot(pcawg_oncogene_sum,
            aes(x = or, y = -log10(p.adj))) +
  geom_point(size = 0.5, alpha = 0.1) +
  geom_point(data = subset(pcawg_oncogene_sum, gene_id %in% hl), color = "red") +
  ggrepel::geom_text_repel(aes(label = gene_id),
                           data = subset(pcawg_oncogene_sum, gene_id %in% hl),
                           max.overlaps = 100) +
  geom_vline(xintercept = tcga_ref, linetype = 2) +
  cowplot::theme_cowplot()  +
  labs(x = "Odds ratio")
  
p2
# MYC 在PCAWG中不显著

ggsave("plots/gene_tendency_TCGA.pdf", p1, width = 7, height = 5)
ggsave("plots/gene_tendency_PCAWG.pdf", p2, width = 7, height = 5)


devtools::load_all("~/proj/gcap/")
devtools::load_all("~/proj/gcaputils//")

# Gene level KM analysis
can_IDs = TCGA_SNP$getGeneSummary()[circular > 2 & gene_id %in% hl_tcga]$gene_id
mat = TCGA_SNP$getGeneSummary(TRUE)
colnames(mat) = gsub("\\.", "-", colnames(mat))
dt_mat = data.table::as.data.table(t(mat), keep.rownames = "sample")
dt_mat = merge(dt_mat, TCGA_SNP$sample_summary[, list(sample, OS.time, OS)], by = "sample")
rv = list()
for (i in can_IDs) {
  dd = dt_mat[, c("sample", i, "OS.time", "OS"), with = FALSE]
  dd = dd[!is.na(dd[[2]])]
  if (nrow(dd) > 5) {
    message("process ", i)
    dd[[2]] = factor(dd[[2]], c("noncircular", "circular"))
    rv[[i]] = ezcox::ezcox(dd, covariates = i, time = "OS.time", status = "OS")
  } else {
    message("skip ", i)
  }
}
rv = rbindlist(rv)
# rv[, fdr := p.adjust(p.value, method = "BH")]
rv2  = rv[p.value < 0.05]

ggpubr::ggdotchart(rv2 %>% dplyr::mutate(
  grp = ifelse(HR > 1, "risk", "protect"),
  x = paste0(Variable, " (", n_contrast, "/", n_ref, ")")
), x = "x", y = "HR",
           color = "grp",                                # Color by groups
           palette = c("#00AFBB",  "#FC4E07"), # Custom color palette
           sorting = "descending",                       # Sort value in descending order
           add = "segments",                             # Add segments from y = 0 to dots
           rotate = TRUE,                                # Rotate vertically
           dot.size = 7,                                 # Large dot size
           label = round(rv2$HR, 2),                        # Add mpg values as dot labels
           font.label = list(color = "white", size = 7, 
                             vjust = 0.5),               # Adjust label parameters
           ggtheme = ggpubr::theme_pubr()                        # ggplot2 theme
) + theme(legend.position = "none") + labs(y = "Hazard ratio", x = NULL) -> p
p
ggsave("plots/oncogene_signif_with_OS.pdf", p, width = 5, height = 6)

rv_tcga = copy(rv)

# pcawg
can_IDs = PCAWG$getGeneSummary()[circular > 2 & gene_id %in% hl_pcawg]$gene_id
mat = PCAWG$getGeneSummary(TRUE)
dt_mat = data.table::as.data.table(t(mat), keep.rownames = "sample")
dt_mat = merge(dt_mat, PCAWG$sample_summary[, list(sample, OS.time, OS)], by = "sample")
rv = list()
for (i in can_IDs) {
  dd = dt_mat[, c("sample", i, "OS.time", "OS"), with = FALSE]
  dd = dd[!is.na(dd[[2]])]
  if (nrow(dd) > 5) {
    message("process ", i)
    dd[[2]] = factor(dd[[2]], c("noncircular", "circular"))
    rv[[i]] = ezcox::ezcox(dd, covariates = i, time = "OS.time", status = "OS")
  } else {
    message("skip ", i)
  }
}
rv = rbindlist(rv)
# rv[, fdr := p.adjust(p.value, method = "BH")]
rv2  = rv[p.value < 0.05 & n_contrast >= 2]

ggpubr::ggdotchart(rv2 %>% dplyr::mutate(
  grp = ifelse(HR > 1, "risk", "protect"),
  x = paste0(Variable, " (", n_contrast, "/", n_ref, ")")
), x = "x", y = "HR",
color = "grp",                                # Color by groups
palette = c("#FC4E07"), # Custom color palette
sorting = "descending",                       # Sort value in descending order
add = "segments",                             # Add segments from y = 0 to dots
rotate = TRUE,                                # Rotate vertically
dot.size = 7,                                 # Large dot size
label = round(rv2$HR, 2),                        # Add mpg values as dot labels
font.label = list(color = "white", size = 7, 
                  vjust = 0.5),               # Adjust label parameters
ggtheme = ggpubr::theme_pubr()                        # ggplot2 theme
) + theme(legend.position = "none") + labs(y = "Hazard ratio", x = NULL) -> p
p
ggsave("plots/oncogene_signif_with_OS_PCAWG.pdf", p, width = 5, height = 3)

rv_pcawg = copy(rv)
openxlsx::write.xlsx(
  list(tcga = rv_tcga, pcawg = rv_pcawg), file = "data/pancan_selected_OS_associated_oncogenes.xlsx"
)

## Circular freq vs Noncircular Freq -----------
# options(IDConverter.datapath = system.file("extdata", package = "IDConverter"))
# 
# data_tcga = copy(TCGA_SNP$data)
# data_tcga[, gene_id := convert_hm_genes(gene_id)]
# data_tcga = data_tcga[!is.na(gene_id)]
# 
# N = nrow(TCGA_SNP$sample_summary)
# data_gene = data_tcga[, list(N = .N, CN = mean(total_cn, na.rm = TRUE)), by = list(gene_id, gene_class)]
# data_gene
# gene_N = dcast(data_gene, gene_id ~ gene_class, value.var = "N", fill = NA)
# gene_N = gene_N[!is.na(noncircular) & !is.na(circular)]
# gene_CN = dcast(data_gene, gene_id ~ gene_class, value.var = "CN", fill = NA)
# gene_CN = gene_CN[!is.na(noncircular) & !is.na(circular)]
# 
# #gene_N = gene_N[circular > quantile(gene_N$circular, 0.5)]
# gene_N = gene_N[circular > 10]
# gene_CN = gene_CN[gene_id %in% gene_N$gene_id]
# gene_dt = merge(gene_N, gene_CN, by = "gene_id")
# colnames(gene_dt) = c("gene_id", "noncircular_freq", "circular_freq", "noncircular_CN", "circular_CN")
# 
# library(ggplot2)
# ggplot(gene_dt, aes(x = noncircular_freq / N, y = circular_freq / N)) +
#   geom_line()
