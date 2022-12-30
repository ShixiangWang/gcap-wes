library(gcap)
library(data.table)
library(IDConverter)
options(IDConverter.datapath = system.file("extdata", package = "IDConverter"))

setwd(file.path(PROJ_DIR, "manuscript"))

TCGA = readRDS("data/TCGA_SNP.rds")
PCAWG = readRDS("data/PCAWG.rds")

info_tcga = TCGA$sample_summary
info_pcawg = PCAWG$sample_summary


# TCGA
SBS = data.table::fread("/data3/wsx_data/MutationalSignatures/TCGA_WES_sigProfiler_SBS_signatures_in_samples.csv")

dt_SBS = info_tcga %>% 
  dplyr::select(sample, class) %>% 
  dplyr::left_join(SBS %>% 
                     dplyr::mutate(sample = substr(`Sample Names`, 1, 15)) %>% 
                     dplyr::select(-c("Cancer Types", "Accuracy", "Sample Names")) %>% 
                     dplyr::distinct(sample, .keep_all = TRUE), by = "sample") %>% 
  dplyr::filter(!is.na(SBS1))

dt_SBS$TMB = rowSums(dt_SBS[, -c(1, 2)])
dt_SBS[, TMB := log2(TMB + 1e-6)]
dt_SBS[, class := factor(class, c("nofocal", "noncircular", "circular"))]

source("../lib/plot.R")

p = plot_clean_stat_box(dt_SBS[, list(class, TMB)], x = class, y = TMB, ylab = "TMB (log2 based)")
ggsave(filename = "plots/TCGA_TMB_by_mutsig.pdf", p, width = 5, height = 4)

dt_SBS2 = data.table::copy(dt_SBS)
dt_SBS2[, MSI := SBS6 + SBS15 + SBS21 + SBS26 + SBS44 + SBS9 + SBS10a + SBS10b + SBS14 + SBS20]
dt_SBS2[, MSI := ifelse(MSI > 0, "Hypermutated signature", "No")]

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
ggsave(filename = "plots/TCGA_hypermut_by_mutsig.pdf", p, width = 4, height = 4)

# PCAWG
SBS = data.table::fread("/data3/wsx_data/MutationalSignatures/PCAWG_sigProfiler_SBS_signatures_in_samples.csv")


dt_SBS = info_pcawg %>% 
  dplyr::select(sample, class) %>% 
  dplyr::left_join(SBS %>% 
                     dplyr::mutate(sample = `Sample Names`) %>% 
                     dplyr::select(-c("Cancer Types", "Accuracy", "Sample Names")) %>% 
                     dplyr::distinct(sample, .keep_all = TRUE), by = "sample") %>% 
  dplyr::filter(!is.na(SBS1))

dt_SBS$TMB = rowSums(dt_SBS[, -c(1, 2)])
dt_SBS[, TMB := log2(TMB + 1e-6)]
dt_SBS[, class := factor(class, c("nofocal", "noncircular", "circular"))]

source("../lib/plot.R")

p = plot_clean_stat_box(dt_SBS[, list(class, TMB)], x = class, y = TMB, ylab = "TMB (log2 based)")
p
ggsave(filename = "plots/PCAWG_TMB_by_mutsig.pdf", p, width = 5, height = 4)

dt_SBS2 = data.table::copy(dt_SBS)
dt_SBS2[, MSI := SBS6 + SBS15 + SBS21 + SBS26 + SBS44 + SBS9 + SBS10a + SBS10b + SBS14 + SBS20]
dt_SBS2[, MSI := ifelse(MSI > 0, "Hypermutated signature", "No")]
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
ggsave(filename = "plots/PCAWG_hypermut_by_mutsig.pdf", p, width = 4, height = 4)
