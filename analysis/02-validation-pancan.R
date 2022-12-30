PROJ_DIR = "~/gcap-analysis/manuscript/"
setwd(PROJ_DIR)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Pancan validation by analysis
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
library(gcap)
library(gcaputils)

# Survival analysis
# TCGA
proj_files = list.files("/data3/wsx_data/tcga_snp_array_gcap_result2", 
                        pattern = "sample_info", full.names = TRUE, all.files = TRUE)
sample_info = purrr::map_df(proj_files, data.table::fread)

amp_files = list.files("/data3/wsx_data/tcga_snp_array_gcap_result2", 
                       pattern = "fCNA_records", full.names = TRUE, all.files = TRUE)
amp_list = purrr::map_df(amp_files, data.table::fread)

TCGA_SNP = fCNA$new(amp_list, sample_info, min_prob = 0.9)

data_snp = dplyr::left_join(
  TCGA_SNP$sample_summary,
  UCSCXenaShiny::tcga_surv,
) %>% 
  dplyr::left_join(UCSCXenaShiny::tcga_clinical %>% 
                     dplyr::select(sample, type))

dput(data_snp[is.na(type), ]$sample)
# check from gdc portal and internet
# 有些是根据临近的barcode推测的
dd = data.table::data.table(
  sample = 
  c("TCGA-08-0373-01", "TCGA-08-0385-01", "TCGA-12-0653-01", "TCGA-12-1601-01", 
    "TCGA-13-1479-01", "TCGA-16-1048-01", "TCGA-28-2501-01", "TCGA-28-2510-01", 
    "TCGA-2G-AALF-01", "TCGA-2G-AALG-01", "TCGA-2G-AALN-01", "TCGA-2G-AALO-01", 
    "TCGA-2G-AALQ-01", "TCGA-2G-AALR-01", "TCGA-2G-AALS-01", "TCGA-2G-AALT-01", 
    "TCGA-2G-AALW-01", "TCGA-2G-AALX-01", "TCGA-2G-AALY-01", "TCGA-2G-AALZ-01", 
    "TCGA-2G-AAM2-01", "TCGA-2G-AAM3-01", "TCGA-2G-AAM4-01", "TCGA-32-2498-01", 
    "TCGA-33-4579-01", "TCGA-36-2539-01", "TCGA-5M-AAT5-01", "TCGA-5M-AATA-01", 
    "TCGA-AN-A0FE-01", "TCGA-AN-A0FG-01", "TCGA-BH-A0B2-01", "TCGA-BR-4186-01", 
    "TCGA-BR-4190-01", "TCGA-BR-4195-01", "TCGA-BR-4196-01", "TCGA-BR-4197-01", 
    "TCGA-BR-4199-01", "TCGA-BR-4200-01", "TCGA-BR-4205-01", "TCGA-BR-4259-01", 
    "TCGA-BR-4260-01", "TCGA-BR-4261-01", "TCGA-BR-4264-01", "TCGA-BR-4265-01", 
    "TCGA-BR-4266-01", "TCGA-BR-4270-01", "TCGA-BR-4271-01", "TCGA-BR-4272-01", 
    "TCGA-BR-4273-01", "TCGA-BR-4274-01", "TCGA-BR-4276-01", "TCGA-BR-4277-01", 
    "TCGA-BR-4278-01", "TCGA-BR-4281-01", "TCGA-BR-4282-01", "TCGA-BR-4283-01", 
    "TCGA-BR-4284-01", "TCGA-BR-4285-01", "TCGA-BR-4286-01", "TCGA-BR-4288-01", 
    "TCGA-BR-4291-01", "TCGA-BR-4298-01", "TCGA-BR-4375-01", "TCGA-BR-4376-01", 
    "TCGA-F4-6857-01", "TCGA-F5-6810-01", "TCGA-O2-A5IC-01", "TCGA-R8-A6YH-01"
  ),
  type = c("GBM", "GBM", "GBM", "GBM", "OV", "GBM", "GBM", "GBM",
           "TGCT", "TGCT", "TGCT", "TGCT", "TGCT", "TGCT", "TGCT", "TGCT",
           "TGCT", "TGCT", "TGCT", "TGCT", "TGCT", "TGCT", "TGCT", "GBM",
           "LUSC", "OV", "COAD", "COAD", "BRCA", "BRCA", "BRCA", "STAD",
           "STAD", "STAD", "STAD", "STAD", "STAD", "STAD", "STAD", "STAD",
           "STAD", "STAD", "STAD", "STAD", "STAD", "STAD", "STAD", "STAD",
           "STAD", "STAD", "STAD", "STAD", "STAD", "STAD", "STAD", "STAD",
           "STAD", "STAD", "STAD", "STAD", "STAD", "STAD", "STAD", "STAD",
           "COAD", "READ", "LUSC", "LGG")
)

data_snp = dplyr::left_join(
  TCGA_SNP$sample_summary,
  UCSCXenaShiny::tcga_surv,
) %>% 
  dplyr::left_join(UCSCXenaShiny::tcga_clinical %>% dplyr::select(sample, type) %>%
                     dplyr::bind_rows(dd)) %>% 
  dplyr::left_join(UCSCXenaShiny::tcga_clinical %>%
                     dplyr::select(
                       sample,
                       age = age_at_initial_pathologic_diagnosis, gender, stage = ajcc_pathologic_tumor_stage))

data_snp = data_snp %>% 
  dplyr::mutate(stage = stringr::str_match(.data$stage, "Stage\\s+(.*?)[ABC]?$")[, 2]) 


sum(is.na(data_snp$type))
saveRDS(data_snp, file = "data/tcga_snp_data.rds")

TCGA_SNP$sample_summary = data_snp
saveRDS(TCGA_SNP, file = "data/TCGA_SNP.rds")

data_snp = readRDS("data/tcga_snp_data.rds")
#undebug(gcap.plotDistribution)
#devtools::load_all("~/proj/gcaputils")
p = gcap.plotDistribution(data_snp[, .(sample, class, by = type)])
p
ggsave("plots/TCGA_fCNA_cancer_types.pdf", p, width = 8, height = 3)

p1 = gcap.plotKMcurve(data_snp[, .(sample, class, OS.time, OS)])
p1
p2 = gcap.plotKMcurve(data_snp[, .(sample, class, PFI.time, PFI)])
p2

pdf("plots/TCGA_SNP_OS_kmplot_for_fCNA_class.pdf", width = 7, height = 7, onefile = FALSE)
print(p1)
dev.off()

pdf("plots/TCGA_SNP_PFI_kmplot_for_fCNA_class.pdf", width = 7, height = 7, onefile = FALSE)
print(p2)
dev.off()

library(survival)
library(forestmodel)
data_snp$class = gcaputils::set_default_factor(data_snp$class)

fit = coxph(Surv(OS.time, OS) ~ class + type, data = data_snp)
p = forestmodel::forest_model(
  fit, covariates = "class", 
  format_options= forestmodel::forest_model_format_options(point_size = 3),
  breaks = c(log(1), log(1.2), log(1.5)))
p
ggplot2::ggsave("plots/TCGA_SNP_OS_forestplot.pdf", p, width = 7, height = 3)

fit = coxph(Surv(PFI.time, PFI) ~ class + type, data = data_snp)
p = forestmodel::forest_model(
  fit, covariates = "class", 
  format_options= forestmodel::forest_model_format_options(point_size = 3),
  breaks = c(log(1), log(1.2), log(1.5)))
p
ggplot2::ggsave("plots/TCGA_SNP_PFS_forestplot.pdf", p, width = 7, height = 3)

summary(coxph(Surv(OS.time, OS) ~ class + type, data = data_snp[class != "nofocal"]))
# p=7.95e-06 HR 7.998e-01
summary(coxph(Surv(PFI.time, PFI) ~ class + type, data = data_snp[class != "nofocal"]))
# p=0.000577 HR 8.490e-01

sum(is.na(data_snp$stage))
sum(is.na(data_snp$age))
sum(is.na(data_snp$gender))

library(ezcox)
t1 = ezcox_group(data_snp, grp_var = "type", covariate = "class",
            controls = c("age", "gender"), time = "OS.time", status = "OS")

p = t1$data$stats %>% 
  dplyr::filter(!is_control & !is.na(HR) & n_contrast > 0) %>% 
  dplyr::mutate(Variable = paste0("(", Group, ") ", contrast_level, " (N=", n_contrast, ")")) %>% 
  dplyr::arrange(contrast_level, desc(HR)) %>% 
  forester(xlim = c(0, 10))

ggplot2::ggsave("plots/TCGA_SNP_OS_type_forest_all.pdf", p, width = 7, height = 12)

plist = list()
for (i in unique(t1$data$stats$Group)) {
  #i = "KICH"
  pf = subset(t1$data$stats, Group == i)
  pval = pf$p.value[pf$contrast_level == "circular"]
  pval2 = pf$p.value[pf$contrast_level == "noncircular"]
  if (isTRUE(pval < 0.05) || isTRUE(pval2 < 0.05)) {
    mr = ifelse(max(pf$upper_95, na.rm = TRUE) > 12, 12, max(pf$upper_95, na.rm = TRUE) + 0.1)
    plist[[i]] = forester(pf, 
                          xlim = c(0, mr),
                          display_cols = c("contrast_level", "HR", "lower_95", "upper_95")) +
      ggplot2::ggtitle(i)
  } else {
    message("skip ", i)
  }
}
length(plist)

p = cowplot::plot_grid(plotlist = plist, nrow = 3)
ggplot2::ggsave("plots/TCGA_SNP_OS_type_forest.pdf", p, width = 12, height = 8)

# add stage
t1 = ezcox_group(data_snp, grp_var = "type", covariate = "class",
                 controls = c("age", "gender", "stage"), time = "OS.time", status = "OS")

p = t1$data$stats %>% 
  dplyr::filter(!is_control & !is.na(HR) & n_contrast > 0) %>% 
  dplyr::mutate(Variable = paste0("(", Group, ") ", contrast_level, " (N=", n_contrast, ")")) %>% 
  dplyr::arrange(contrast_level, desc(HR)) %>% 
  forester(xlim = c(0, 10))

ggplot2::ggsave("plots/TCGA_SNP_OS_type_forest_all_stage.pdf", p, width = 7, height = 10)

plist = list()
for (i in unique(t1$data$stats$Group)) {
  #i = "KICH"
  pf = subset(t1$data$stats, Group == i)
  pval = pf$p.value[pf$contrast_level == "circular"]
  pval2 = pf$p.value[pf$contrast_level == "noncircular"]
  if (isTRUE(pval < 0.05) || isTRUE(pval2 < 0.05)) {
    mr = ifelse(max(pf$upper_95, na.rm = TRUE) > 12, 12, max(pf$upper_95, na.rm = TRUE) + 0.1)
    plist[[i]] = forester(pf, 
                          xlim = c(0, mr),
                          display_cols = c("contrast_level", "HR", "lower_95", "upper_95")) +
      ggplot2::ggtitle(i)
  } else {
    message("skip ", i)
  }
}
length(plist)

p = cowplot::plot_grid(plotlist = plist, nrow = 1)
ggplot2::ggsave("plots/TCGA_SNP_OS_type_forest2.pdf", p, width = 16, height = 4)

## PFI
t1 = ezcox_group(data_snp, grp_var = "type", covariate = "class",
                 controls = c("age", "gender"), time = "PFI.time", status = "PFI")

plist = list()
for (i in unique(t1$data$stats$Group)) {
  #i = "KICH"
  pf = subset(t1$data$stats, Group == i)
  pval = pf$p.value[pf$contrast_level == "circular"]
  pval2 = pf$p.value[pf$contrast_level == "noncircular"]
  if (isTRUE(pval < 0.05) || isTRUE(pval2 < 0.05)) {
    mr = ifelse(max(pf$upper_95, na.rm = TRUE) > 12, 12, max(pf$upper_95, na.rm = TRUE) + 0.1)
    plist[[i]] = forester(pf, 
                          xlim = c(0, mr),
                          display_cols = c("contrast_level", "HR", "lower_95", "upper_95")) +
      ggplot2::ggtitle(i)
  } else {
    message("skip ", i)
  }
}

length(plist)
p = cowplot::plot_grid(plotlist = plist, nrow = 2)
ggplot2::ggsave("plots/TCGA_SNP_PFS_type_forest.pdf", p, width = 16, height = 5)

# add stage
t1 = ezcox_group(data_snp, grp_var = "type", covariate = "class",
                 controls = c("age", "gender", "stage"), time = "PFI.time", status = "PFI")
plist = list()
for (i in unique(t1$data$stats$Group)) {
  #i = "KICH"
  pf = subset(t1$data$stats, Group == i)
  pval = pf$p.value[pf$contrast_level == "circular"]
  pval2 = pf$p.value[pf$contrast_level == "noncircular"]
  if (isTRUE(pval < 0.05) || isTRUE(pval2 < 0.05)) {
    mr = ifelse(max(pf$upper_95, na.rm = TRUE) > 12, 12, max(pf$upper_95, na.rm = TRUE) + 0.1)
    plist[[i]] = forester(pf, 
                          xlim = c(0, mr),
                          display_cols = c("contrast_level", "HR", "lower_95", "upper_95")) +
      ggplot2::ggtitle(i)
  } else {
    message("skip ", i)
  }
}
length(plist)

p = cowplot::plot_grid(plotlist = plist, nrow = 1)
ggplot2::ggsave("plots/TCGA_SNP_PFS_type_forest2.pdf", p, width = 12, height = 4)


# Survival analysis
# PCAWG
proj_files = list.files("/data3/wsx_data/pcawg_gcap_result2", 
                        pattern = "sample_info", full.names = TRUE, all.files = TRUE)
sample_info = purrr::map_df(proj_files, data.table::fread)

amp_files = list.files("/data3/wsx_data/pcawg_gcap_result2/", pattern = "fCNA_records", full.names = TRUE, all.files = TRUE)
amp_list = purrr::map_df(amp_files[file.info(amp_files)$size > 100], data.table::fread)

sample_info2 = merge(readRDS("../pancan-analysis/data/pcawg_sample_info.rds")[
  , .(sample, donor_id, gender = donor_sex, age = donor_age_at_diagnosis, 
      donor_survival_time, donor_vital_status, 
      dcc_project_code, cancer_type)], sample_info, by = "sample")
sample_info2[, `:=`(OS.time = donor_survival_time, OS = as.integer(donor_vital_status == "deceased"))]
sample_info2[, type := stringr::str_split(cancer_type, "\\-", simplify = T)[, 1]]

PCAWG = fCNA$new(amp_list, sample_info2, min_prob = 0.9)
saveRDS(PCAWG, file = "data/PCAWG.rds")

PCAWG = readRDS("data/PCAWG.rds")
data_pcawg = PCAWG$sample_summary

table(data_pcawg$class)

data_pcawg[, class := set_default_factor(class)]
data_pcawg

devtools::load_all("~/proj/gcaputils")
p = gcap.plotDistribution(data_pcawg[, .(sample, class, by = cancer_type)])
#p = gcap.plotDistribution(data_pcawg[, .(sample, class, by = type)])
p
ggsave("plots/PCAWG_fCNA_cancer_types.pdf", p, width = 8, height = 3)

p1 = gcap.plotKMcurve(data_pcawg[, .(sample, class, OS.time, OS)])
p1
pdf("plots/PCAWG_OS_kmplot_for_fCNA_class.pdf", width = 7, height = 7, onefile = FALSE)
print(p1)
dev.off()

fit = coxph(Surv(OS.time, OS) ~ class + type, data = data_pcawg)
p = forestmodel::forest_model(
  fit, covariates = "class", 
  format_options= forestmodel::forest_model_format_options(point_size = 3),
  breaks = c(log(1), log(1.2), log(1.5), log(2)))
p
ggplot2::ggsave("plots/PCAWG_OS_forestplot.pdf", p, width = 7, height = 3)

summary(coxph(Surv(OS.time, OS) ~ class + type, data = data_pcawg[class != "nofocal"]))
# p=0.048295 HR 0.83

t1 = ezcox_group(data_pcawg, grp_var = "type", covariate = "class",
                 controls = c("age", "gender"), time = "OS.time", status = "OS")

p = t1$data$stats %>% 
  dplyr::filter(!is_control & !is.na(HR) & n_contrast > 0) %>% 
  dplyr::mutate(Variable = paste0("(", Group, ") ", contrast_level, " (N=", n_contrast, ")")) %>% 
  dplyr::arrange(contrast_level, desc(HR)) %>% 
  forester(xlim = c(0, 10))

ggplot2::ggsave("plots/PCAWG_OS_type_forest_all.pdf", p, width = 7, height = 8)

plist = list()
for (i in unique(t1$data$stats$Group)) {
  #i = "KICH"
  pf = subset(t1$data$stats, Group == i)
  pval = pf$p.value[pf$contrast_level == "circular"]
  pval2 = pf$p.value[pf$contrast_level == "noncircular"]
  if (isTRUE(pval < 0.05) || isTRUE(pval2 < 0.05)) {
    mr = ifelse(max(pf$upper_95, na.rm = TRUE) > 12, 12, max(pf$upper_95, na.rm = TRUE) + 0.1)
    plist[[i]] = forester(pf, 
                          xlim = c(0, mr),
                          display_cols = c("contrast_level", "HR", "lower_95", "upper_95")) +
      ggplot2::ggtitle(i)
  } else {
    message("skip ", i)
  }
}

length(plist)
p = cowplot::plot_grid(plotlist = plist, nrow = 1)
ggplot2::ggsave("plots/PCAWG_OS_type_forest.pdf", p, width = 16, height = 3)

# Mutational signature
# TCGA
library(dplyr)
SBS = data.table::fread("/data3/wsx_data/MutationalSignatures/TCGA_WES_sigProfiler_SBS_signatures_in_samples.csv")
dt_SBS = data_snp %>% 
  dplyr::select(sample, class) %>% 
  dplyr::left_join(SBS %>% 
                     dplyr::mutate(sample = substr(`Sample Names`, 1, 15)) %>% 
                     dplyr::select(-c("Cancer Types", "Accuracy", "Sample Names")) %>% 
                     dplyr::distinct(sample, .keep_all = TRUE), by = "sample") %>% 
  dplyr::filter(!is.na(SBS1))

CNS = readxl::read_excel("/data3/wsx_data/Nature2022_CNS/Nature2022_CNS.xlsx", sheet = 3)
info = readxl::read_excel("/data3/wsx_data/Nature2022_CNS/Nature2022_ASCAT.xlsx", sheet = 3) %>% data.table()


dt_CNS = data_snp %>% 
  dplyr::select(sample, class) %>% 
  dplyr::left_join(
    dplyr::left_join(info %>% dplyr::select(name, barcodeTumour), 
                     CNS, by = c("name" = "Sample")) %>% 
      dplyr::mutate(sample = substr(barcodeTumour, 1, 15)) %>% 
      dplyr::select(-c("barcodeTumour", "name")) %>% 
      dplyr::mutate_at(dplyr::vars(-sample), as.numeric), by = "sample")
which(is.na(dt_CNS))


call_grp_analysis = function(dt, ref_group = NA) {
  library(sigminer)
  ge = group_enrichment(dt,
                        grp_vars = "class",
                        enrich_vars = colnames(dt)[-c(1, 2)],
                        co_method = "wilcox.test",
                        ref_group = ref_group)
  if (is.na(ref_group)) ge$grp1 = set_default_factor(ge$grp1)
  ge$enrich_var = factor(ge$enrich_var, colnames(dt)[-c(1, 2)])
  ge$fdr = ifelse(ge$fdr == 0, ge$fdr + .Machine$double.xmin, ge$fdr)
  p = show_group_enrichment(ge, cut_p_value = TRUE, cluster_row = FALSE, return_list = TRUE)
  p = p[[1]] + labs(x = NULL, y = NULL) + ggpubr::rotate_x_text() #+ coord_flip()
  return(list(data = ge, plot = p))
}

rv_SBS = call_grp_analysis(dt_SBS)
rv_SBS$plot

rv_CNS = call_grp_analysis(dt_CNS)
rv_CNS$plot

dt_CNS2 = data_snp[, c("sample", "class", paste0("CN", 1:19))]
rv_CNS2 = call_grp_analysis(dt_CNS2)
rv_CNS2$plot

ggsave("plots/tcga_fCNA_SBS.pdf", rv_SBS$plot,
       width = 22, height = 2.5, device = cairo_pdf)
ggsave("plots/tcga_fCNA_CNS_rel_from_paper.pdf", rv_CNS$plot,
       width = 10, height = 2.5, device = cairo_pdf)
ggsave("plots/tcga_fCNA_CNS_abs_from_sigminer.pdf", rv_CNS2$plot,
       width = 10, height = 2.5, device = cairo_pdf)

# Shorten the SBS by removing non-sigficant data
dt_SBS2 = dt_SBS[, c("sample", "class",
                     as.character(unique(rv_SBS$data[fdr < 0.05]$enrich_var))), with = FALSE]
rv_SBS2 = call_grp_analysis(dt_SBS2)
rv_SBS2$plot

ggsave("plots/tcga_fCNA_SBS_v2.pdf", rv_SBS2$plot,
       width = 15, height = 2.5, device = cairo_pdf)

if (file.exists("data/MutationalSignatureComparison_TCGA.xlsx")) file.remove("data/MutationalSignatureComparison_TCGA.xlsx")
openxlsx::write.xlsx(
  list(
    SBS = rv_SBS$data,
    CNS_rel = rv_CNS$data,
    CNS_abs = rv_CNS2$data
  ), file = "data/MutationalSignatureComparison_TCGA.xlsx"
)

source("../lib/plot.R")

dt_APOBEC = dt_SBS[, .(sample, class, apobec = SBS2+SBS13)]
dt_APOBEC

ggpubr::ggviolin(dt_APOBEC, x = "class", y = "apobec") +
  ggpubr::stat_compare_means(method = "wilcox.test")

p = ggstatsplot::ggbetweenstats(
  dt_APOBEC, # %>% dplyr::mutate(apobec = ifelse(apobec > 1000, 1000, apobec)),
  x = class,
  y = apobec,
  bf.message = FALSE, type = "np", centrality.type = "p"
)

p = p + cowplot::theme_cowplot() + theme(legend.position = "none") +
  ylab("APOBEC associated mutations")
p

#plot_clean_stat_box(dt_APOBEC, x = class, y = apobec, ylab = "APOBEC associated mutations")
ggsave("plots/APOBEC_signature_activity_TCGA.pdf", p, width = 8, height = 5)

p = ggstatsplot::ggbetweenstats(
  dt_CNS,
  x = class,
  y = CN8,
  bf.message = FALSE, type = "np", centrality.type = "p"
)
p = p + cowplot::theme_cowplot() + theme(legend.position = "none") +
  ylab("CN8 contribution")
p

ggsave("plots/CN8_rel_activity_TCGA.pdf", p, width = 8, height = 5)

p = ggstatsplot::ggbetweenstats(
  dt_CNS2,
  x = class,
  y = CN8,
  bf.message = FALSE, type = "np", centrality.type = "p"
)
p = p + cowplot::theme_cowplot() + theme(legend.position = "none") +
  ylab("CN8 contribution (abs)")
p

ggsave("plots/CN8_abs_activity_TCGA.pdf", p, width = 8, height = 5)

library(sigminer)
library(IDConverter)
options(IDConverter.datapath = system.file("extdata", package = "IDConverter"))

apobec_pcawg = data.table::fread("/data3/wsx_data/MutationalSignatures/MAF_Aug31_2016_sorted_A3A_A3B_comparePlus.sp")
apobec_pcawg = merge(apobec_pcawg, data_pcawg[, list(sample, class)], by.x = "Tumor_Sample_Barcode", by.y = "sample")
source("../lib/plot.R")
p = plot_clean_stat_box(apobec_pcawg, class, APOBECtCa_enrich, ylab = "APOBECtCa_enrich")
ggsave("plots/PCAWG_APOBECtCa_enrich.pdf", width = 5, height = 4)

#dt_APOBEC %>% dplyr::group_by(class) %>% dplyr::summarise(m = mean(apobec))

# Mutational signature
# PCAWG

library(dplyr)
SBS = data.table::fread("/data3/wsx_data/MutationalSignatures/PCAWG_sigProfiler_SBS_signatures_in_samples.csv")
DBS = data.table::fread("/data3/wsx_data/MutationalSignatures/PCAWG_sigProfiler_DBS_signatures_in_samples.csv")
ID = data.table::fread("/data3/wsx_data/MutationalSignatures/PCAWG_SigProfiler_ID_signatures_in_samples.csv")

dt_CNS = data_pcawg[, c("sample", "class", paste0("CN", 1:19))]

dt_SBS = data_pcawg %>% 
  dplyr::select(sample, class) %>% 
  dplyr::left_join(SBS %>% 
                     dplyr::mutate(sample = `Sample Names`) %>% 
                     dplyr::select(-c("Cancer Types", "Accuracy", "Sample Names")) %>% 
                     dplyr::distinct(sample, .keep_all = TRUE), by = "sample") %>% 
  dplyr::filter(!is.na(SBS1))

dt_DBS = data_pcawg %>% 
  dplyr::select(sample, class) %>% 
  dplyr::left_join(DBS %>% 
                     dplyr::mutate(sample = `Sample Names`) %>% 
                     dplyr::select(-c("Cancer Types", "Accuracy", "Sample Names")) %>% 
                     dplyr::distinct(sample, .keep_all = TRUE), by = "sample")

dt_ID = data_pcawg %>% 
  dplyr::select(sample, class) %>% 
  dplyr::left_join(ID %>% 
                     dplyr::mutate(sample = `Sample Names`) %>% 
                     dplyr::select(-c("Cancer Types", "Accuracy", "Sample Names")) %>% 
                     dplyr::distinct(sample, .keep_all = TRUE), by = "sample")


rv_SBS = call_grp_analysis(dt_SBS)
rv_SBS$plot

rv_DBS = call_grp_analysis(dt_DBS)
rv_DBS$plot

rv_ID = call_grp_analysis(dt_ID)
rv_ID$plot

ggsave("plots/pcawg_fCNA_SBS.pdf", rv_SBS$plot,
       width = 22, height = 2.5, device = cairo_pdf)
ggsave("plots/pcawg_fCNA_DBS.pdf", rv_DBS$plot,
       width = 6, height = 2.5, device = cairo_pdf)
ggsave("plots/pcawg_fCNA_ID.pdf", rv_ID$plot,
       width = 9, height = 2.5, device = cairo_pdf)

# Shorten the SBS by removing non-sigficant data
dt_SBS2 = dt_SBS[,  c("sample", "class",
                      as.character(unique(rv_SBS$data[fdr < 0.05]$enrich_var))), with = FALSE]
rv_SBS2 = call_grp_analysis(dt_SBS2)
rv_SBS2$plot

ggsave("plots/pcawg_fCNA_SBS_v2.pdf", rv_SBS2$plot,
       width = 10, height = 2.5, device = cairo_pdf)


rv_CNS = call_grp_analysis(dt_CNS)
rv_CNS$plot

ggsave("plots/pcawg_fCNA_CNS.pdf", rv_CNS$plot,
       width = 10, height = 2.5, device = cairo_pdf)

dt_APOBEC = dt_SBS[, .(sample, class, apobec = SBS2+SBS13)]
dt_APOBEC

p = ggstatsplot::ggbetweenstats(
  dt_APOBEC,
  x = class,
  y = apobec,
  bf.message = FALSE, type = "np", centrality.type = "p"
)

p = p + cowplot::theme_cowplot() + theme(legend.position = "none") +
  ylab("APOBEC associated mutations")
p
ggsave("plots/APOBEC_signature_activity_PCAWG.pdf", p, width = 8, height = 5)

dt_CN8 = dt_CNS[, .(sample, class, CN8)]
dt_CN8$CN8 = dt_CN8$CN8 / rowSums(dt_CNS[, -c(1, 2)])
dt_CN8$class = set_default_factor(dt_CN8$class)

# Follow the operation from Nature paper
dt_CN8[, CN8 := ifelse(CN8 <= 0.05, 0, CN8)]

p = ggstatsplot::ggbetweenstats(
  dt_CN8,
  x = class,
  y = CN8,
  bf.message = FALSE, type = "np", centrality.type = "p"
)
p = p + cowplot::theme_cowplot() + theme(legend.position = "none") +
  ylab("CN8 contribution")
p

ggsave("plots/CN8_rel_activity_PCAWG.pdf", p, width = 8, height = 5)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Model performance
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
perf_gene = readRDS("../pancan-analysis/data/gene_perf_of_model_on_TCGA_SNP_data.rds")
perf_gene = perf_gene[model %in% c("xgb11", "CN8", "CN7_plus_CN8")]
perf_gene$model[1] = "XGB"
colnames(perf_gene)[2:3] = c("auPRC", "auROC")
perf_gene = melt(perf_gene, id.vars = "model", variable.name = "measure", value.name = "score")

p = perf_gene %>% 
  ggplot(aes(x = model, y = score)) + 
  geom_col(fill = "steelblue") +
  facet_wrap(~measure, nrow = 1) +
  geom_text(aes(label = round(score, 3)), nudge_y = 0.05) +
  cowplot::theme_cowplot() +
  xlab(NULL) + ylab("Performance score")
p = p + ggpubr::rotate_x_text(30)
p
ggsave("plots/TCGA_snp_model_gene_performance.pdf", p, width = 7, height = 4)

perf = readRDS("../pancan-analysis/data/sample_perf_of_model_on_TCGA_SNP_data.rds")
perf = perf[model %in% c("xgb11", "CN8", "CN7_plus_CN8")]
perf$model[1] = "XGB"
colnames(perf)[2:3] = c("auPRC", "auROC")
perf = melt(perf, id.vars = "model", variable.name = "measure", value.name = "score")

p = perf %>% 
  ggplot(aes(x = model, y = score)) + 
  geom_col(fill = "steelblue") +
  facet_wrap(~measure, nrow = 1) +
  geom_text(aes(label = round(score, 3)), nudge_y = 0.05) +
  cowplot::theme_cowplot() +
  xlab(NULL) + ylab("Performance score")
p = p + ggpubr::rotate_x_text(30)
p
ggsave("plots/TCGA_snp_model_sample_performance.pdf", p, width = 7, height = 4)

perf2 = rbind(cbind(perf_gene, type = "gene"),
              cbind(perf, type = "sample"))[measure %in% c("auPRC", "auROC")]

perf2

p = perf2 %>% 
  ggplot(aes(x = measure, y = score)) + 
  geom_col(aes(fill = model), position = position_dodge2()) +
  facet_wrap(~type, nrow = 1) +
  geom_text(aes(y = score + 0.05, label = round(score, 3), group = model),
            position = position_dodge(width = 0.9), size = 3) +
  cowplot::theme_cowplot() + theme(legend.position = "top") +
  xlab(NULL) + ylab("Performance score") + scale_fill_brewer()
#p = p + ggpubr::rotate_x_text(30)
p
ggsave("plots/TCGA_snp_model_all_performance.pdf", p, width = 7, height = 4)


