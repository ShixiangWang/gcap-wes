library(gcap)
library(data.table)

setwd(file.path(PROJ_DIR, "pancan-analysis"))


# TCGA --------------------------------------------------------------------

proj_files = list.files("/data3/wsx_data/tcga_snp_array_gcap_result", pattern = "sample_info", full.names = TRUE, all.files = TRUE)
sample_info = purrr::map_df(proj_files, data.table::fread)

amp_files = list.files("/data3/wsx_data/tcga_snp_array_gcap_result", pattern = "fCNA_records", full.names = TRUE, all.files = TRUE)
amp_list = purrr::map_df(amp_files, data.table::fread)

TCGA_SNP = fCNA$new(amp_list, sample_info)
saveRDS(TCGA_SNP, file = "data/TCGA_SNP.rds")

data_snp = TCGA_SNP$sample_summary

# Combine circular and possibly_circular
data_snp[, class := fcase(class %in% c("circular", "possibly_circular"), "circular",
                           class == "noncircular", "noncircular",
                           default = "nofocal")]
data_snp[, class := factor(class, c("nofocal", "noncircular", "circular"))]


table(data_snp$class)
#https://gdc.cancer.gov/resources-tcga-users/tcga-code-tables/sample-type-codes
table(type = substr(data_snp$sample, 14, 15), class = data_snp$class)

data_snp = dplyr::left_join(
  data_snp[as.integer(substr(sample, 14, 15)) < 10, .(sample = substr(sample, 1, 15), class)],
  UCSCXenaShiny::tcga_surv,
) %>% 
  dplyr::left_join(UCSCXenaShiny::tcga_clinical %>% dplyr::select(sample, type))


library(survival)

fit = survfit(Surv(OS.time, OS) ~ class, data = data_snp)
p = survminer::ggsurvplot(fit, pval = TRUE, data = data_snp,
                          palette = c("grey", "blue", "red"),
                          risk.table = TRUE, ylab = "OS",
                          legend.labs = levels(data_snp$class))
p

pdf("plots/TCGA_SNP_OS_kmplot_for_fCNA_class.pdf", width = 7, height = 7, onefile = FALSE)
print(p)
dev.off()

fit = survfit(Surv(PFI.time, PFI) ~ class, data = data_snp)
p = survminer::ggsurvplot(fit, pval = TRUE, data = data_snp,
                          palette = c("grey", "blue", "red"),
                          risk.table = TRUE, ylab = "PFS",
                          legend.labs = levels(data_snp$class) #, censor.size = 2, size = 0.2
                          )
p

pdf("plots/TCGA_SNP_PFS_kmplot_for_fCNA_class.pdf", width = 7, height = 7, onefile = FALSE)
print(p)
dev.off()

summary(coxph(Surv(OS.time, OS) ~ class + strata(type), data = data_snp))
library(see)
p = plot(parameters::parameters(coxph(Surv(OS.time, OS) ~ class + strata(type), data = data_snp), exponentiate = T),
     show_labels  = TRUE)
p
ggsave("plots/TCGA_SNP_OS_simplified_forestplot.pdf", plot = p, width = 5, height = 3)

p = plot(parameters::parameters(coxph(Surv(PFI.time, PFI) ~ class + strata(type), data = data_snp), exponentiate = T),
         show_labels  = TRUE)
p
ggsave("plots/TCGA_SNP_PFS_simplified_forestplot.pdf", plot = p, width = 5, height = 3)


fit = coxph(Surv(OS.time, OS) ~ class + type, data = data_snp)
p = forestmodel::forest_model(fit, covariates = "class", format_options= forestmodel::forest_model_format_options(point_size = 3))
p
ggplot2::ggsave("plots/TCGA_SNP_OS_forestplot.pdf", p, width = 7, height = 3)


fit = coxph(Surv(PFI.time, PFI) ~ class + type, data = data_snp)
p = forestmodel::forest_model(fit, covariates = "class", format_options= forestmodel::forest_model_format_options(point_size = 3))
p
ggplot2::ggsave("plots/TCGA_SNP_PFS_forestplot.pdf", p, width = 7, height = 3)


summary(coxph(Surv(OS.time, OS) ~ class + strata(type), data = data_snp[class != "nofocal"]))
summary(coxph(Surv(PFI.time, PFI) ~ class + strata(type), data = data_snp[class != "nofocal"]))

# PCAWG -------------------------------------------------------------------

#proj_files = list.files("/data3/wsx_data/pcawg_gcap_result/", pattern = "sample_info", full.names = TRUE, all.files = TRUE)
#sample_info = purrr::map_df(proj_files, data.table::fread)

amp_files = list.files("/data3/wsx_data/pcawg_gcap_result/", pattern = "fCNA_records", full.names = TRUE, all.files = TRUE)
amp_list = purrr::map_df(amp_files[file.info(amp_files)$size > 100], data.table::fread)

PCAWG = fCNA$new(amp_list, readRDS("data/pcawg_sample_info.rds"))
saveRDS(PCAWG, file = "data/PCAWG.rds")
data_pcawg = PCAWG$sample_summary

table(data_pcawg$class)

data_pcawg[, `:=`(OS.time = donor_survival_time, OS = as.integer(donor_vital_status == "deceased"))]
data_pcawg[, class := data.table::fcase(class %in% c("circular", "possibly_circular"), "circular",
                            class == "noncircular", "noncircular",
                            default = "nofocal")]
data_pcawg[, class := factor(class, c("nofocal", "noncircular", "circular"))]
data_pcawg[, type := stringr::str_split(cancer_type, "\\-", simplify = T)[, 1]]

# 仔细看看了NG 文章的图，只有按下面手动截尾方式处理才能得到相似的区间和最后不低到0的生存概率
#data_pcawg[, OS := ifelse(OS.time >= 365 * 5, 0, OS)]
#data_pcawg[, OS.time := ifelse(OS.time >= 365 * 5, 365 * 5, OS.time)]


library(survival)

fit = survfit(Surv(OS.time, OS) ~ class, data = data_pcawg)
p = survminer::ggsurvplot(fit, pval = TRUE, data = data_pcawg,
                          palette = c("grey", "blue", "red"),
                          risk.table = TRUE, ylab = "OS",
                          legend.labs = levels(data_pcawg$class))
p

pdf("plots/PCAWG_OS_kmplot_for_fCNA_class.pdf", width = 7, height = 7, onefile = FALSE)
print(p)
dev.off()

summary(coxph(Surv(OS.time, OS) ~ class + strata(type), data = data_pcawg))

library(see)
p = plot(parameters::parameters(coxph(Surv(OS.time, OS) ~ class + strata(type), data = data_pcawg), exponentiate = T),
         show_labels  = TRUE)
p
ggsave("plots/PCAWG_OS_simplified_forestplot.pdf", plot = p, width = 5, height = 3)

fit = coxph(Surv(OS.time, OS) ~ class + type, data = data_pcawg)
p = forestmodel::forest_model(fit, covariates = "class", format_options= forestmodel::forest_model_format_options(point_size = 3))
p
ggplot2::ggsave("plots/PCAWG_OS_forestplot.pdf", p, width = 7, height = 3)


