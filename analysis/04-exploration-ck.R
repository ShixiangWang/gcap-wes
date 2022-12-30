PROJ_DIR = "~/gcap-analysis/manuscript/"
setwd(PROJ_DIR)

library(data.table)
library(dplyr)
library(gcap)

case_result = data.table::fread("/data3/wsx_data/ck_gcap_result/CK1000_sample_info.csv")
clinfo = data.table::fread("../SYSUCC-CRC/data/SYSUCC_CRC_Clinical_withDFS_Mutannotated.txt")
clinfo[, Tumor_Sample_Barcode := stringr::str_remove(Tumor_Sample_Barcode, "_NT")]
clinfo = unique(clinfo[, .(Tumor_Sample_Barcode, CancerType, gender,
                           age = age_at_diagnosis,
                           smoking,
                           primary_tumor_location,
                           metastasis = metastasis_at_diagnosis,
                           pathological_stage,
                           TMB,
                           MSI_Status = MSIstatus, family_history,
                           Cluster,
                           OS, OS_Status = Osstatus, DFS = DFS_Day, DFS_Status)])
clinfo[, OS_Status := ifelse(OS_Status == "1", 1L, 0L)]


data_case = dplyr::left_join(clinfo, case_result, by = c("Tumor_Sample_Barcode" = "sample"))
colnames(data_case)[1] = "sample"

table(data_case$class)

CRC = fCNA$new(fcna = data.table::fread("/data3/wsx_data/ck_gcap_result/CK1000_fCNA_records.csv"), 
               pdata = data_case)

table(CRC$sample_summary$class)

saveRDS(CRC, file = "data/SYSUCC_CRC.rds")

CRC = readRDS("data/SYSUCC_CRC.rds")

# Sample analysis -------------------------------------------------------
library(survival)
library(gcaputils)

p = gcap.plotKMcurve(CRC, surv_data = c("OS", "OS_Status"))
p

pdf("plots/SYSUCC_CRC_OS_kmplot_for_fCNA_class.pdf", width = 7, height = 6, onefile = FALSE)
print(p)
dev.off()


p = gcap.plotKMcurve(CRC, surv_data = c("DFS", "DFS_Status"))
p

pdf("plots/SYSUCC_CRC_DFS_kmplot_for_fCNA_class.pdf", width = 7, height = 6, onefile = FALSE)
print(p)
dev.off()

data_case = CRC$sample_summary
data_case$class = set_default_factor(data_case$class)
data_case$Cluster = factor(data_case$Cluster, c("GS", "CIN-LR", "CIN-HR", "Hypermutated"))
data_case[, pathological_stage := ifelse(pathological_stage == "unknown", NA, pathological_stage)]

str(data_case$class)
table(data_case$class)

data.table::setnames(data_case, c("class", "Cluster"), c("fCNA_cluster", "SYSUCC_cluster"))

cor(data_case$ploidy, data_case$cna_burden, use = "pairwise")
# [1] 0.9070007

fit = coxph(Surv(OS, OS_Status) ~ fCNA_cluster + SYSUCC_cluster + 
              MSI_Status + TMB +
              CancerType + family_history +
              smoking + primary_tumor_location +
              gender + age + pathological_stage,
            data = data_case)

p = forestmodel::forest_model(fit, format_options= forestmodel::forest_model_format_options(point_size = 3))
p
ggplot2::ggsave("plots/SYSUCC_CRC_OS_forestplot_for_multivariables.pdf", p, width = 9, height = 8)

fit = coxph(Surv(DFS, DFS_Status) ~ fCNA_cluster + SYSUCC_cluster + 
              MSI_Status + TMB +
              CancerType + family_history +
              smoking + primary_tumor_location +
              gender + age + pathological_stage,
            data = data_case)
p = forestmodel::forest_model(fit, format_options= forestmodel::forest_model_format_options(point_size = 3))
p
ggplot2::ggsave("plots/SYSUCC_CRC_DFS_forestplot_for_multivariables.pdf", p, width = 9, height = 8)

table(data_case$SYSUCC_cluster, data_case$fCNA_cluster)

source("../lib/plot.R")
data_case$cnv = rowSums(data_case[, paste0("CN", 1:19)])

p = plot_clean_stat_box(data_case %>% mutate(cnv = log2(cnv + 1e-6)), x = fCNA_cluster, y = cnv,
                    ylab = "CN segments (log2 based)")
p
ggsave("plots/SYSUCC_CRC_CNV_events_by_class.pdf", p, width = 5, height = 4)

p = plot_clean_stat_box(data_case, x = fCNA_cluster, y = cna_burden,
                        ylab = "CNA burden")

p
ggsave("plots/SYSUCC_CRC_CNA_burden_by_class.pdf", p, width = 5, height = 4)

p = plot_clean_stat_box(data_case %>% mutate(TMB = log2(TMB + 1e-6)), x = fCNA_cluster, y = TMB,
                        ylab = "TMB (log2 based)")

p
ggsave("plots/SYSUCC_CRC_TMB_by_class.pdf", p, width = 5, height = 4)

p = plot_clean_stat_box(data_case, x = fCNA_cluster, y = pLOH,
                        ylab = "Genome percentage with LOH")

p
ggsave("plots/SYSUCC_CRC_pLOH_by_class.pdf", p, width = 5, height = 4)

p = plot_clean_stat_box(data_case, x = fCNA_cluster, y = AScore,
                        ylab = "Aneuploidy score")

p
ggsave("plots/SYSUCC_CRC_Aneuploidy_score_by_class.pdf", p, width = 5, height = 4)

p = plot_clean_stat_box(data_case, x = fCNA_cluster, y = purity,
                        ylab = "Tumor purity")

p
ggsave("plots/SYSUCC_CRC_purity_by_class.pdf", p, width = 5, height = 4)

p = plot_clean_stat_box(data_case, x = fCNA_cluster, y = ploidy,
                        ylab = "Tumor ploidy")

p
ggsave("plots/SYSUCC_CRC_ploidy_by_class.pdf", p, width = 5, height = 4)


grp_ca = get_group_comparison(data_case, col_group = "fCNA_cluster",
                              cols_to_compare = c("CancerType",
                                                  "gender",
                                                  "metastasis",
                                                  "pathological_stage",
                                                  "smoking",
                                                  "primary_tumor_location",
                                                  "MSI_Status",
                                                  "SYSUCC_cluster"))
plist = show_group_comparison(grp_ca, ca_p_threshold = 0.001, 
                              legend_position_ca = "top", legend_title_ca = "class",
                              xlab = NA,
                              set_ca_sig_yaxis = T, text_angle_x = 0, text_hjust_x = 0.5)
plist$ca_comb

plotlist = lapply(plist$ca, function(x) x + 
                    ggplot2::scale_fill_manual(values = c("grey", "blue", "red"),
                                               name = "class"))

p = cowplot::plot_grid(plotlist = plotlist[c(1:3, 5:7)], 
  align = "hv", ncol = 3)
p
p2 = cowplot::plot_grid(plotlist = plotlist[c(4, 8)], 
                       align = "hv", ncol = 1)
p2
p3 = cowplot::plot_grid(p, p2, ncol = 2, rel_widths = c(5, 3))
p3
ggsave("plots/SYSUCC_CRC_CAs_by_class.pdf", p3, width = 13, height = 6)

library(ComplexHeatmap)
library(circlize)
dt = data_case[!is.na(CN1)]
dt = dt[!is.na(SYSUCC_cluster)]
dt[, combined := paste0(fCNA_cluster, "/", SYSUCC_cluster)]
# > table(dt$combined)
# 
# circular/CIN-HR          circular/CIN-LR              circular/GS    circular/Hypermutated           nofocal/CIN-HR 
# 49                       56                       57                        2                      108 
# nofocal/CIN-LR               nofocal/GS     nofocal/Hypermutated       noncircular/CIN-HR       noncircular/CIN-LR 
# 71                      346                       69                       63                      125 
# noncircular/GS noncircular/Hypermutated 
# 56                        1 

dt[, combined := ifelse(SYSUCC_cluster == "Hypermutated", "Hypermutated", combined)]
table(dt$combined)

library(survminer)

dt$combined = factor(dt$combined)
dt$combined = relevel(dt$combined, ref = "Hypermutated")
dt$c = dt$combined
dt$c = factor(dt$c, levels = c("Hypermutated",
                               "nofocal/CIN-LR",
                               "nofocal/GS",
                               "noncircular/CIN-LR", "noncircular/GS", "noncircular/CIN-HR",
                               "circular/GS", "circular/CIN-LR", "nofocal/CIN-HR",
                               "circular/CIN-HR"))
p = forestmodel::forest_model(coxph(Surv(OS, OS_Status) ~ c, data = dt),
                              format_options =
                                forestmodel::forest_model_format_options(point_size = 3))
p
ggplot2::ggsave("plots/SYSUCC_CRC_OS_combined_cluster.pdf", p, width = 7, height = 5)


dt[, cluster2 := fcase(
  combined == "Hypermutated", "HM",
  combined %in% "nofocal/CIN-LR", "CIN-Mild",
  combined %in% "nofocal/GS", "CN-Quiet",
  combined %in% c("noncircular/GS", "noncircular/CIN-LR"), "Non-Circ",
  combined %in% c("nofocal/CIN-HR", "noncircular/CIN-HR", "circular/GS", "circular/CIN-LR"), "CIN-HR|Circ",
  combined == "circular/CIN-HR", "CIN-HR&Circ"
)]

new_lvls = c("HM", "CIN-Mild", "CN-Quiet", "Non-Circ", "CIN-HR|Circ", "CIN-HR&Circ")
dt$cluster2 = factor(dt$cluster2, levels = new_lvls)

m1 = coxph(Surv(OS, OS_Status) ~ SYSUCC_cluster, data = dt)
m2 = coxph(Surv(OS, OS_Status) ~ SYSUCC_cluster + fCNA_cluster, data = dt)
m3 = coxph(Surv(OS, OS_Status) ~ fCNA_cluster, data = dt)
m1$concordance["concordance"]
# 0.5744264 
m2$concordance["concordance"]
# 0.6034302 
m3$concordance["concordance"]
# 0.5601656 
anova(m1, m2)
# 0.01702 *
anova(m3, m2)
# 0.0003117 ***

m1 = coxph(Surv(DFS, DFS_Status) ~ SYSUCC_cluster, data = dt)
m2 = coxph(Surv(DFS, DFS_Status) ~ SYSUCC_cluster + fCNA_cluster, data = dt)
m3 = coxph(Surv(DFS, DFS_Status) ~ fCNA_cluster, data = dt)
m1$concordance["concordance"]
# 0.5852873
m2$concordance["concordance"]
# 0.6150337
m3$concordance["concordance"]
# 0.5793017
anova(m1, m2)
# 0.0006736 ***
anova(m3, m2)
# 0.001495 **


md = data.frame(
  variable = c("m1", "m2", "m1", "m2"),
  value = c(0.5744264, 0.6034302, 0.5852873, 0.6150337),
  g2 = c("OS", "OS", "DFS", "DFS")
)
md$g2 = factor(md$g2, levels = c("OS", "DFS"))


ggplot(md, 
       aes(x = variable, y = value, fill = variable)) +
  geom_bar(stat = "identity", width = 0.5) +
  facet_grid(. ~ g2) +
  cowplot::theme_cowplot() +
  ylim(0, 0.7) + labs(x = NULL, y = "C index") +
  theme(legend.position = "none") -> p
p
ggsave("plots/SYSUCC_CK_c_index.pdf", p, width = 4, height = 3)

p = forestmodel::forest_model(coxph(Surv(OS, OS_Status) ~ cluster2, data = dt),
                              format_options= forestmodel::forest_model_format_options(point_size = 3))
p
ggplot2::ggsave("plots/SYSUCC_CRC_OS_cluster2.pdf", p, width = 7, height = 4)

new_pal_cols = c("grey", "black", "cyan", "blue", "magenta", "red")

fit = surv_fit(Surv(OS, OS_Status) ~ cluster2, data = dt)
cls_lvls = sub("cluster2=", "", names(fit$strata))
pdf("plots/SYSUCC_CRC_OS_KM_cluster2.pdf", width = 6, height = 6, onefile = FALSE)
ggsurvplot(fit, 
           palette = new_pal_cols,
           risk.table = TRUE, risk.table.height = 0.35,
           legend.labs = cls_lvls, size = 0.5, censor.size = 2,
           pval = TRUE)
dev.off()

fit = surv_fit(Surv(DFS, DFS_Status) ~ cluster2, data = dt)
cls_lvls = sub("cluster2=", "", names(fit$strata))

pdf("plots/SYSUCC_CRC_DFS_KM_cluster2.pdf", width = 6, height = 6, onefile = FALSE)
ggsurvplot(fit, 
           palette = new_pal_cols,
           risk.table = TRUE, risk.table.height = 0.35,
           legend.labs = cls_lvls, size = 0.5, censor.size = 2,
           pval = TRUE)
dev.off()

source("../lib/plot.R")
dt$cnv = rowSums(dt[, paste0("CN", 1:19)])

colors = new_pal_cols
p = plot_clean_stat_box(dt %>% mutate(cnv = log2(cnv + 1e-6)), x = cluster2, y = cnv,
                        ylab = "CN segments (log2 based)", colors)
p
ggsave("plots/SYSUCC_CRC_CNV_events_by_cluster2.pdf", p, width = 7, height = 6)

p = plot_clean_stat_box(dt, x = cluster2, y = cna_burden,
                        ylab = "CNA burden", colors)

p
ggsave("plots/SYSUCC_CRC_CNA_burden_by_cluster2.pdf", p, width = 6, height = 6)

p = plot_clean_stat_box(dt %>% mutate(TMB = log2(TMB + 1e-6)), x = cluster2, y = TMB,
                        ylab = "TMB (log2 based)", colors)

p
ggsave("plots/SYSUCC_CRC_TMB_by_cluster2.pdf", p, width = 7, height = 6)

p = plot_clean_stat_box(dt, x = cluster2, y = pLOH,
                        ylab = "Genome percentage with LOH", colors)

p
ggsave("plots/SYSUCC_CRC_pLOH_by_cluster2.pdf", p, width = 7, height = 6)

p = plot_clean_stat_box(dt, x = cluster2, y = AScore,
                        ylab = "Aneuploidy score", colors)

p
ggsave("plots/SYSUCC_CRC_Aneuploidy_score_by_cluster2.pdf", p, width = 7, height = 6)

p = plot_clean_stat_box(dt, x = cluster2, y = purity,
                        ylab = "Tumor purity", colors)

p
ggsave("plots/SYSUCC_CRC_purity_by_cluster2.pdf", p, width = 7, height = 5)

p = plot_clean_stat_box(dt, x = cluster2, y = ploidy,
                        ylab = "Tumor ploidy", colors)

p
ggsave("plots/SYSUCC_CRC_ploidy_by_cluster2.pdf", p, width = 7, height = 6)


grp_ca = get_group_comparison(dt, col_group = "cluster2",
                              cols_to_compare = c("gender",
                                                  "metastasis",
                                                  "pathological_stage",
                                                  "primary_tumor_location"))
plist = show_group_comparison(grp_ca, ca_p_threshold = 0.001, 
                              legend_position_ca = "top", legend_title_ca = "class",
                              xlab = NA,
                              set_ca_sig_yaxis = T, text_angle_x = 0, text_hjust_x = 0.5)
plist$ca_comb

plotlist = lapply(plist$ca, function(x) x + 
                    ggplot2::scale_fill_manual(values = colors,
                                               name = "class"))

p = cowplot::plot_grid(plotlist = plotlist, 
                       align = "hv", ncol = 2)
p

ggsave("plots/SYSUCC_CRC_CAs_by_cluster2.pdf", p, width = 10, height = 6)


saveRDS(dt[, list(sample, cluster2)], file = "data/SYSUCC_CRC_sample_cluster2.rds")

# Mutational signature enrichment
sbs = readxl::read_excel("data/SYSUCC_CRC_MutSig.xlsx", 1)
p = readxl::read_excel("data/SYSUCC_CRC_MutSig.xlsx", 4)
kept = p %>% dplyr::filter(N >= 10) %>% dplyr::pull(sig)
library(dplyr)
sbs = sbs %>% 
  filter(sig %in% kept, ! sig %in% paste0("SBS", c("27", "43", 45:60)) ) %>%
  mutate(sample = sub("_NT", "", sample)) %>%
  group_by(sample, sig) %>% 
  summarise(expo = median(expo, na.rm = TRUE)) %>% 
  tidyr::pivot_wider(id_cols = "sample", names_from = "sig", values_from = "expo", values_fill = 0) %>%
  data.table::as.data.table()
sbs = merge(dt2[, list(sample, cluster2)], sbs, by = "sample")
sbs

ms = merge(sbs, dt[, c("sample", paste0("CN", 1:19)), with = FALSE])

ms_grp = group_enrichment(ms, grp_vars = "cluster2",
                          enrich_vars = setdiff(colnames(ms), c("sample", "cluster2")),
                          co_method = "wilcox.test")

ms_grp[, Foldchange := measure_observed]
library(ggplot2)

plot_enrich = function(df, .g = "HM") {
  ggplot(df[grp1 == .g], aes(x = Foldchange, y = -log10(fdr))) +
    geom_point(size = 0.5) +
    geom_vline(xintercept = 1, linetype = 2) +
    geom_hline(yintercept = -log10(0.05), linetype = 2) +
    geom_point(
      color = "red", 
      data = subset(df[grp1 == .g], Foldchange > 1 & fdr < 0.05)
    ) +
    ggrepel::geom_text_repel(
      aes(label = enrich_var), 
      data = subset(df[grp1 == .g], Foldchange > 1 & fdr < 0.05), max.overlaps = 100,
      size = 2
    ) +
    cowplot::theme_cowplot()
}

p1 = plot_enrich(ms_grp) +
  scale_x_log10(breaks = c(1, 10, 100, 1000)) + labs(title = "HM")
p1
p2 = plot_enrich(ms_grp, .g = "CN-Quiet") + labs(title = "CN-Quiet")
p3 = plot_enrich(ms_grp, .g = "CIN-Mild") + labs(title = "CIN-Mild")
p4 = plot_enrich(ms_grp, .g = "Non-Circ") + labs(title = "Non-Circ")
p5 = plot_enrich(ms_grp, .g = "CIN-HR|Circ") + labs(title = "CIN-HR|Circ")
p6 = plot_enrich(ms_grp, .g = "CIN-HR&Circ") + labs(title = "CIN-HR&Circ")

p = cowplot::plot_grid(
  p1, p2, p3, p4, p5, p6, align = "hv", nrow = 2
)
ggsave("plots/SYSUCC_cluster2_grp_enrichment_ms.pdf", plot = p, width = 14, height = 6)

show_order = function(.g) {
  dd = ms_grp[grp1 == .g & fdr < 0.05 & Foldchange > 1]
  dd[order(fdr, -Foldchange), list(enrich_var, Foldchange, fdr)]
}
show_order("HM")$enrich_var
show_order("CIN-Mild")$enrich_var
show_order("CN-Quiet")$enrich_var
show_order("Non-Circ")$enrich_var
show_order("CIN-HR|Circ")$enrich_var
show_order("CIN-HR&Circ")$enrich_var

openxlsx::write.xlsx(ms_grp, file = "data/SYSUCC_CRC_cluster2_MutSig_enrichment.xlsx")

library(maftools)
maf = data.table::fread("../SYSUCC-CRC/data/SYSUCC_CRC_Mutation_withSilent.maf")
maf[, Tumor_Sample_Barcode := sub("_NT", "", Tumor_Sample_Barcode)]
maf = maf[Tumor_Sample_Barcode %in% dt$sample]
maf = read.maf(maf,
               clinicalData = dt[, list(Tumor_Sample_Barcode = sample, cluster2)])

maf_HM = subsetMaf(maf, clinQuery = paste0('cluster2 %in% "', "HM", '"'))
maf_Quiet = subsetMaf(maf, clinQuery = paste0('cluster2 %in% "', "CN-Quiet", '"'))
maf_Others = subsetMaf(maf, clinQuery = paste0('!cluster2 %in% c("HM", "CN-Quiet")'))

OncogenicPathways(maf_HM, )
OncogenicPathways(maf_Quiet)
OncogenicPathways(maf_Others)

PlotOncogenicPathways(maf_HM, "TP53",removeNonMutated = FALSE)
PlotOncogenicPathways(maf_Quiet, "TP53", removeNonMutated = FALSE)
PlotOncogenicPathways(maf_Others, "TP53", removeNonMutated = FALSE)

onco_hm = oncodrive(maf_HM)
pdf("plots/SYSUCC_CRC_maf_HM_oncoplot.pdf", width = 11, height = 8)
oncoplot(maf_HM, genes = onco_hm[fdr < 0.05]$Hugo_Symbol, logColBar = TRUE, removeNonMutated = FALSE)
dev.off()

onco_quiet = oncodrive(maf_Quiet)
pdf("plots/SYSUCC_CRC_maf_quiet_oncoplot.pdf", width = 11, height = 4)
oncoplot(maf_Quiet, genes = onco_quiet[fdr < 0.05]$Hugo_Symbol, logColBar = TRUE, removeNonMutated = FALSE)
dev.off()


onco_others = oncodrive(maf_Others)
pdf("plots/SYSUCC_CRC_maf_others_oncoplot.pdf", width = 11, height = 14)
oncoplot(maf_Others, genes = onco_others[fdr < 0.05]$Hugo_Symbol, logColBar = TRUE, removeNonMutated = FALSE)
dev.off()

save(onco_hm, onco_quiet, onco_others, file = "data/SYSUCC_CRC_oncodrive.RData")

HM_vs_others = maftools::mafCompare(maf_HM, maf_Others)
Quiet_vs_others = maftools::mafCompare(maf_Quiet, maf_Others)

Quiet_vs_others$results[adjPval < 0.05]

# Code below is outdated
mat = t(as.matrix(data.table::setDF(dt2[, c(-1, -2)], rownames = dt2$sample)))
col_fun = colorRamp2(c(0, 0.4), c("white", "red"))

ht = Heatmap(mat, col = col_fun, cluster_columns = FALSE,
        show_column_names = FALSE, column_split = dt2$cluster2,
        heatmap_legend_param = list(title = "Contribution",
                                    color_bar = "discrete", col = "red"),
        top_annotation = HeatmapAnnotation(
          cluster = dt$cluster2,
            col = list(
              cluster = c("HM" = "grey", "L-Risk"="green", "M-Risk"="orange", "H-Risk"="red", "E-Risk"="purple")
            )
        ))
pdf("plots/SYSUCC_cluster2_CN_signature.pdf", width = 9, height = 5)
draw(ht, merge_legends = TRUE)
dev.off()

# Remove CN1, CN9 (diploidy) before clustering
mat = t(as.matrix(data.table::setDF(dt2[, c(-1, -2)], rownames = dt2$sample)))
# mat = mat[!rownames(mat) %in% c("CN1", "CN2", "CN3", "CN18", "CN19"), ]
mat = mat[!rownames(mat) %in% c("CN1"), ]
col_fun = colorRamp2(c(0, 0.2), c("white", "red"))

ht = Heatmap(mat, col = col_fun, cluster_columns = TRUE, cluster_rows = FALSE,
             show_column_names = FALSE, column_split = dt2$cluster2,
             heatmap_legend_param = list(title = "Contribution",
                                         color_bar = "discrete", col = "red"),
             top_annotation = HeatmapAnnotation(
               cluster = dt2$cluster2,
               col = list(
                 cluster = c("HM" = "grey", "L-Risk"="green", "M-Risk"="orange", "H-Risk"="red", "E-Risk"="purple")
               )
             ))
pdf("plots/SYSUCC_cluster2_CN_signature_cluster2.pdf", width = 9, height = 5)
draw(ht, merge_legends = TRUE)
dev.off()

mat = t(as.matrix(data.table::setDF(sbs[, c(-1, -2)], rownames = sbs$sample)))
col_fun = colorRamp2(c(0, 0.4), c("white", "red"))

ht = Heatmap(mat, col = col_fun, cluster_columns = FALSE,
             show_column_names = FALSE, column_split = sbs$cluster2,
             heatmap_legend_param = list(title = "Contribution",
                                         color_bar = "discrete", col = "red"),
             top_annotation = HeatmapAnnotation(
               cluster = sbs$cluster2,
               col = list(
                 cluster = c("HM" = "grey", "L-Risk"="green", "M-Risk"="orange", "H-Risk"="red", "E-Risk"="purple")
               )
             ))
pdf("plots/SYSUCC_cluster2_SBS_signature.pdf", width = 9, height = 6)
draw(ht, merge_legends = TRUE)
dev.off()

dt3 = merge(dt2, sbs[, -c("cluster2")], by = "sample")
mat = t(as.matrix(data.table::setDF(dt3[, c(-1, -2)], rownames = dt3$sample)))
col_fun = colorRamp2(c(0, 0.2), c("white", "red"))

ht = Heatmap(mat, col = col_fun, cluster_columns = TRUE,
             show_column_names = FALSE, column_split = dt3$cluster2,
             heatmap_legend_param = list(title = "Contribution",
                                         color_bar = "discrete", col = "red"),
             top_annotation = HeatmapAnnotation(
               cluster = dt3$cluster2,
               col = list(
                 cluster = c("HM" = "grey", "L-Risk"="green", "M-Risk"="orange", "H-Risk"="red", "E-Risk"="purple")
               )
             ))
draw(ht)

dir.create("plots/SYSUCC_MS")
for (i in colnames(ms)[-c(1, 2)]) {
  p = plot_clean_stat_box(ms, x = cluster2, y = !!i,
                          ylab = i, colors = colors)
  p
  ggsave(sprintf("plots/SYSUCC_MS/%s_by_cluster2.pdf", i), p, width = 8, height = 6)
}

file.remove("data/SYSUCC_CRC_cluster2_MutSig.xlsx")
openxlsx::write.xlsx(
  ms, file = "data/SYSUCC_CRC_cluster2_MutSig.xlsx"
)

# Gene analysis -----------------------------------------------------------

CRC = readRDS("data/SYSUCC_CRC.rds")
View(CRC$sample_summary)
table(CRC$sample_summary$pathological_stage)
CRC$sample_summary[, pathological_stage := ifelse(pathological_stage == "unknown", NA_character_, pathological_stage)]


CRC$sample_summary %>% 
  tidyr::unite("subtypes", 
               primary_tumor_location, pathological_stage,
               na.rm = TRUE, remove = FALSE) %>% 
  dplyr::mutate(
    subtypes = ifelse(is.na(pathological_stage), NA_character_, subtypes)
  ) -> t

CRC$sample_summary = data.table::as.data.table(t)
table(t$subtypes)

summary(coxph(Surv(OS, OS_Status) ~ CancerType, data = t))
summary(coxph(Surv(OS, OS_Status) ~ gender, data = t))
summary(coxph(Surv(OS, OS_Status) ~ primary_tumor_location, data = t))
summary(coxph(Surv(OS, OS_Status) ~ CancerType + primary_tumor_location + gender + pathological_stage, data = t))

source("../lib/hl.R")

devtools::load_all("~/proj/gcaputils/")
pdf("plots/SYSUCC_CRC_genome_heatmap.pdf", width = 5, height = 8)
gcap.plotGenomeHeatmap(CRC, group = "CancerType", bycol = TRUE, sort = FALSE,
                       highlights = get_highlights(CRC, group = "CancerType",
                                                   gene_freq = 2, grp_freq = 2), fontsize = 6)
dev.off()

pdf("plots/SYSUCC_CRC_genome_heatmap_oncogene.pdf", width = 4.6, height = 8)
gcap.plotGenomeHeatmap(CRC, group = "CancerType", bycol = TRUE, sort = FALSE,
                       highlights = get_highlights_onco(CRC, group = "CancerType",
                                                   gene_freq = 2, grp_freq = 2), fontsize = 6)
dev.off()

hl_genes = highlights = get_highlights_onco(CRC, group = "CancerType",
                                            gene_freq = 1, grp_freq = 1, return_simple = TRUE)

?gcap.plotCircos
CRC2 = readRDS("data/SYSUCC_CRC.rds")
CRC2$convertGeneID()
undebug(gcap.plotCircos)
pdf("plots/SYSUCC_CRC_circos.pdf", width = 10, height = 10)
gcap.plotCircos(CRC2, genome_build = "hg19",
                highlights = get_highlights_onco(CRC, group = "CancerType",
                                                 gene_freq = 1, grp_freq = 1, return_simple = TRUE))
dev.off()

devtools::load_all("~/proj/gcap/")
devtools::load_all("~/proj/gcaputils//")
hl_genes
CRC2$getCytobandSummary(unique = TRUE)
head(CRC2$getCytobandSummary(unique = TRUE), 20)

can_IDs = CRC2$getCytobandSummary(unique = TRUE)[circular >= 5]$band
rvplot = list()
for (i in can_IDs) {
  rvplot[[i]] = gcap.plotKMcurve(
    CRC2, c("OS", "OS_Status"),
    ID = i,
    mat = CRC2$getCytobandSummary(return_mat = TRUE),
    focus = "all"
  )
  rvplot[[i]]$plot = rvplot[[i]]$plot + ggtitle(i)
}

pdf("plots/SYSUCC_CRC_band_kmplot_list.pdf", width = 28, height = 6)
survminer::arrange_ggsurvplots(rvplot, ncol = 5)
dev.off()

# CRC = readRDS("data/CRC1000.rds")
# 
# CRC$convertGeneID(genome_build = "hg19")
# gene_info = CRC$gene_summary[!is.na(gene_id)][order(circular, Total, decreasing = TRUE)]
# sum(gene_info$circular > 2)
# na.omit(unique(CRC$data[amplicon_type == "circular" & total_cn > 40]$gene_id))
# 
# gcap.plotProfile(CRC, top_n = 100, show_column_names = FALSE)
# 
# #pdf("plots/CRC_fCNA_profile.pdf", width = 10, height = 6)
# # 这种图像不适合大样本和CN相关数据
# gcap.plotProfile(CRC, genes = na.omit(CRC$gene_summary[(circular + possibly_circular) > 2]$gene_id),
#                  show_column_names = FALSE, remove_empty_columns = FALSE,
#                  only_circular = FALSE, merge_circular = TRUE, show_row_names = FALSE)
# #dev.off()
# 
# gcap.plotDistribution(CRC, x = gene_info[circular > 2]$gene_id)
# 
# gcap.plotCircos(CRC, genome_build = "hg19")
# 
# genes_summary = CRC$data[
#   , .(cn = mean(total_cn[amplicon_type %in% c("circular", "possibly_circular")], na.rm = TRUE),
#       prob = mean(prob[amplicon_type %in% c("circular", "possibly_circular")], na.rm = TRUE),
#       N = sum(amplicon_type %in% c("circular", "possibly_circular"), na.rm = TRUE)
#   ), by = gene_id][!is.na(gene_id) & N > 0]
# 
# pdf("plots/CRC_circos.pdf", width = 10, height = 10)
# gcap.plotCircos(CRC, genome_build = "hg19", 
#                 highlight_genes = genes_summary[(cn > 20 & N > 1) | cn > 40][, .(gene_id, label = N)])
# dev.off()
# 
# p = gcap.dotplot(CRC, max.overlaps = 100, color = "red", 
#                  filter = cn > 80 | (N > 1 & cn > 40) | (N > 4 & cn > 20),
#                  size = 2, by = "gene_id")
# p
# ggsave("plots/CRC_dotplot.pdf", p, width = 7, height = 6)
# 
# p = gcap.dotplot(CRC, max.overlaps = 100, color = "red", 
#                  filter = cn > 80 | (N > 1 & cn > 40) | (N > 4 & cn > 20) | N > 50, size = 2, by = "band")
# p
# ggsave("plots/CRC_dotplot_by_cytoband.pdf", p, width = 7, height = 6)
# 
# p = gcap.dotplot(CRC, max.overlaps = 100, color = "red", filter = N > 0, size = 2, by = "chr")
# p
# ggsave("plots/CRC_dotplot_by_chr.pdf", p, width = 7, height = 6)
# 
# #CRC$data[startsWith(band, "chr17") & amplicon_type == "circular"]
# 
# circ_cls = clusterGPosition(gene_info[circular > 0][, .(gene_id)], no_annotation = TRUE, genome_build = "hg19", simplify = FALSE)
# 
# # circular vs noncircular for same gene cluster
# zz = gcap.plotForest(CRC, response_data = c("OS", "OS_Status"),
#                      x_is_gene = TRUE,
#                      x = circ_cls[, .(cluster = paste0("cluster", cluster), gene_id)],
#                      gene_focus = "vs")
# 
# zz$model$plot_forest(ref_line = 1, xlim = c(0, 10))
# circ_cls[cluster == 11]
# 
# # circular vs (noncircular and nofocal)
# zz2 = gcap.plotForest(CRC, response_data = c("OS", "OS_Status"),
#                       x_is_gene = TRUE,
#                       x = circ_cls[, .(cluster = paste0("cluster", cluster), gene_id)],
#                       gene_focus = "circular", optimize_model = TRUE)
# zz2$data
# zz2$optmodel
# zz2$model$plot_forest(ref_line = 1, xlim = c(0, 10), p = 0.05)
# # zz2$model$plot_forest(ref_line = 1, xlim = c(0, 10), vars = names(zz2$optmodel$xlevels))
# 
# # 
# # circ_cls[cluster == 5]
# # CRC$gene_summary[gene_id %in% circ_cls[cluster == 5]$gene_id]
# # 
# # circ_cls2 = data.table(
# #   cluster = "genes",
# #   gene_id = circ_cls[, .(cluster = paste0("cluster", cluster), gene_id)][cluster %in% names(zz2$optmodel$xlevels)]$gene_id
# # )
# # zz3 = gcap.plotForest(CRC, response_data = c("OS", "OS_Status"),
# #                       x_is_gene = TRUE,
# #                       x = circ_cls2,
# #                       gene_focus = "circular")
# # zz3$data

# 
# # ORA ---------------------------------------------------------------------
# library(clusterProfiler)
# CRC = readRDS(file = "data/CRC1000.rds")
# 
# ORA_hallmark = gcap.enrich(CRC, target = "circular")
# ORA_reactome = gcap.enrich(CRC, target = "circular", category = "C2", subcategory = "CP:REACTOME")
# ORA_KEGG = gcap.enrich(CRC, target = "circular", category = "C2", subcategory = "CP:KEGG")
# ORA_BP = gcap.enrich(CRC, target = "circular", category = "C5", subcategory = "GO:BP")
# ORA_MF = gcap.enrich(CRC, target = "circular", category = "C5", subcategory = "GO:MF")
# ORA_CC = gcap.enrich(CRC, target = "circular", category = "C5", subcategory = "GO:CC")
# ORA_TFT = gcap.enrich(CRC, target = "circular", category = "C3", subcategory = "TFT:GTRD")
# ORA_immune = gcap.enrich(CRC, target = "circular",  category = "C7", subcategory = "IMMUNESIGDB")
# 
# dotplot(ORA_hallmark, showCategory = 20, label_format = 50)
# dotplot(ORA_reactome, label_format = 50, showCategory = 20)
# dotplot(ORA_KEGG, showCategory = 20, label_format = 50)
# dotplot(ORA_BP, showCategory = 20, label_format = 50)
# dotplot(ORA_MF, showCategory = 20, label_format = 50)
# dotplot(ORA_CC, showCategory = 20, label_format = 50)
# dotplot(ORA_TFT, showCategory = 20, label_format = 50)
# dotplot(ORA_immune, showCategory = 20, label_format = 50)
# 
# ORA_BP@result = ORA_BP@result %>% 
#   dplyr::mutate(
#     Description = stringr::str_remove(Description, "GOBP_") %>% stringr::str_to_title() %>% stringr::str_replace_all("_", " ")
#   )
# 
# p = dotplot(ORA_BP, showCategory = 40, label_format = 80)
# p
# ggsave("plots/CRC_GO_BP.pdf", p, width = 10, height = 10)
# 
# ORA_reactome@result = ORA_reactome@result %>% 
#   dplyr::mutate(
#     Description = stringr::str_remove(Description, "REACTOME_") %>% stringr::str_to_title() %>% stringr::str_replace_all("_", " ")
#   )
# 
# p = dotplot(ORA_reactome, label_format = 50, showCategory = 20)
# ggsave("plots/CRC_reactome.pdf", p, width = 10, height = 7)
# 
# # > IDConverter::convert_hm_genes(unlist(stringr::str_split(subset(ORA_BP@result, ID == "GOBP_T_CELL_ACTIVATION_INVOLVED_IN_IMMUNE_RESPONSE")$geneID, "/")))
# # Downloading https://zenodo.org/record/6342397/files/human_hg38_gene_info.rds to /tmp/Rtmp7vM668/human_hg38_gene_info.rds
# # trying URL 'https://zenodo.org/record/6342397/files/human_hg38_gene_info.rds'
# # Content type 'application/octet-stream' length 1052986 bytes (1.0 MB)
# # ==================================================
# #   downloaded 1.0 MB
# # 
# # [1] "SMAD7"   "PTGER4"  "HLA-DRA" "HLA-DMB" "CRACR2A" "HMGB1"   "IFNA21"  "IFNA5"   "IFNA16"  "IFNB1"   "IFNW1"   "IFNA10"  "IFNA7"  
# # [14] "IFNA14"  "IFNA17"  "IFNA4"   "RPS6"    "CD74"    "RARA"   

# Gene detection for FISH ----------------

gene_result = data.table::fread("/data3/wsx_data/ck_gcap_result/CK1000_fCNA_records.csv")
library(IDConverter)
options(IDConverter.datapath = system.file("extdata", package = "IDConverter"))
gene_result$gene_name = IDConverter::convert_hm_genes(gene_result$gene_id, genome_build = "hg19")
gene_result

available_samples = c(
  "CRC1030", "CRC1023",
  "CRC1000", "CRC1119",
  "CRC1055", "CRC1001",
  "CRC1022", "CRC1029",
  "CRC1081", "CRC1021",
  "CRC1039", "CRC1044",
  "CRC1077", "CRC1057",
  "CRC1003", "CRC1026",
  "CRC1012", "CRC1002",
  "CRC1106", "CRC1047",
  "CRC1073", "CRC1071",
  "CRC1045", "CRC1048",
  "CRC1017", "CRC1064",
  "CRC1020", "CRC1033",
  "CRC1085", "CRC1123",
  "CRC1122", "CRC1089",
  "CRC1109", "CRC1107"
)

sum(available_samples %in% gene_result$sample)
# 19

keep = gene_result[sample %in% available_samples & gene_class != "nofocal"]
length(unique(keep$sample))
# 19 samples with amplicon detected
# 11 samples with circular amplicon
keep = keep[order(-total_cn, -prob)]
keep_band = keep[gene_class == "circular", list(gene_name = paste(unique(gene_name), collapse = ","),
                        circ_cn_mean = mean(total_cn),
                        prob_mean = mean(prob)), list(sample, band)]
keep_band = keep_band[order(sample, -prob_mean)]

data.table::fwrite(
  keep[, .(sample, band, gene_id, gene_name, total_cn, minor_cn, ploidy, prob)],
  file = "data/SYSUCC_CRC_gene_list_for_available_samples.csv"
)

data.table::fwrite(
  keep_band,
  file = "data/SYSUCC_CRC_cytoband_circular_gene_list_for_available_samples.csv"
)



# Mutational signatures ---------------------------------------------------

library(sigminer)
CRC = readRDS("data/SYSUCC_CRC.rds")
CRC_maf = read_maf("../SYSUCC-CRC/data/SYSUCC_CRC_Mutation_withSilent.maf", verbose = TRUE)
maf_mat = sig_tally(CRC_maf)

# expo_list = sig_fit_bootstrap_batch(t(maf_mat$nmf_matrix),
#                                     sig_db = "SBS", sig_index = "ALL",
#                                     n = 100, use_parallel = 20)
# 
# expo_list$expo
# expo = expo_list$expo[type != "optimal", .(act = median(exposure, na.rm = TRUE)), by = .(sample, sig)]
# expo
# 
# expo[, sample := stringr::str_remove(sample, "_NT")]
# expo
# expo_wide = dcast(expo, sample ~ sig, value.var = "act", fill = NA_real_)
# sig_drop = names(which(sapply(expo_wide[, -1], function(x) max(x, na.rm = TRUE) < 10)))
# 
# expo_wide = expo_wide[, !colnames(expo_wide) %in% sig_drop, with = FALSE]

#
expo_wide = sig_fit(t(maf_mat$nmf_matrix), 
                    sig_index = gsub("SBS", "", rownames(get_sig_db("SBS")$aetiology)[1:65]),
                    sig_db = "SBS",
                    return_class = "data.table")
expo_wide[, sample := stringr::str_remove(sample, "_NT")]
#colnames(expo_wide) = sub("COSMIC_", "SBS", colnames(expo_wide))
expo_wide
#

data = merge(expo_wide, CRC$sample_summary[, .(sample, class)], by = "sample")
data.table::setcolorder(data, c("sample", "class"))

call_grp_analysis = function(dt, ref_group = NA) {
  library(sigminer)
  library(gcaputils)
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

rv_SBS = call_grp_analysis(data)
rv_SBS$plot

ggsave("plots/SYSUCC_CRC_fCNA_SBS.pdf", rv_SBS$plot,
       width = 22, height = 2.5, device = cairo_pdf)

rv_CNS = call_grp_analysis(CRC$sample_summary[, c("sample", "class", paste0("CN", 1:19)), with = FALSE])
rv_CNS$plot

ggsave("plots/SYSUCC_CRC_fCNA_CNS.pdf", rv_CNS$plot,
       width = 8, height = 2.5, device = cairo_pdf)


rv_SBS = call_grp_analysis(data, ref_group = "noncircular")
rv_SBS$plot

ggsave("plots/SYSUCC_CRC_fCNA_SBS_noncircular_ref.pdf", rv_SBS$plot,
       width = 12, height = 2.5, device = cairo_pdf)

rv_CNS = call_grp_analysis(CRC$sample_summary[, c("sample", "class", paste0("CN", 1:19)), with = FALSE], ref_group = "noncircular")
rv_CNS$plot

ggsave("plots/SYSUCC_CRC_fCNA_CNS_noncircular_ref.pdf", rv_CNS$plot,
       width = 8, height = 2.5, device = cairo_pdf)
