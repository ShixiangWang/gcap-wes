library(UCSCXenaShiny)
library(dplyr)
library(readxl)
library(survival)

# Reproduce NG Cox analysis with TCGA data (part of paper data) --------------------------------
sample_tb <- readxl::read_excel("preprocessing-and-EDA/data/ICGC-ecDNA.xlsx", sheet = 3) # The supplementary table

surv_df <- sample_tb %>% 
  dplyr::filter(startsWith(sample_barcode, "TCGA"), tumor_or_normal == "tumor") %>% 
  dplyr::select(sample_barcode, patient_barcode, lineage, sample_classification) %>% 
  dplyr::left_join(tcga_clinical %>%
                     dplyr::select(sample, type), by = c(sample_barcode = "sample")) %>%
  dplyr::left_join(tcga_surv, by = c(sample_barcode = "sample")) %>% 
  dplyr::mutate(sample_classification = factor(
    sample_classification, 
    c("No-fSCNA", "Linear", "Heavily-rearranged", "BFB", "Circular")))

# From paper author:
# When we performed the survival analysis (including cox), 
# we filtered out some samples with high os times (5 yr regardless of their os status).

surv_df = surv_df %>% 
  dplyr::filter(!is.na(OS.time), OS.time <= 365 * 5)

coxph(Surv(OS.time, OS) ~ sample_classification, surv_df)  # This is consistent with paper

# Cox analysis with cancer type controlled
coxph(Surv(OS.time, OS) ~ sample_classification + strata(type), surv_df)
# Use the type label from manuscript
coxph(Surv(OS.time, OS) ~ sample_classification + strata(lineage), surv_df)

# The survival analysis is consistent in only TCGA patients with OS <5yrs.
fit = coxph(Surv(OS.time, OS) ~ sample_classification + strata(type), surv_df)
forestmodel::forest_model(fit)

# What abount COAD
surv_coad = surv_df %>% dplyr::filter(lineage == "Colorectal")
fit = coxph(Surv(OS.time, OS) ~ sample_classification, surv_coad)
summary(fit)
forestmodel::forest_model(fit)

# 单个瘤种中不一定能够观测到相同的趋势

# What about PFI?
surv_df2 <- sample_tb %>% 
  dplyr::filter(startsWith(sample_barcode, "TCGA"), tumor_or_normal == "tumor") %>% 
  dplyr::select(sample_barcode, patient_barcode, lineage, sample_classification) %>% 
  dplyr::left_join(tcga_clinical %>%
                     dplyr::select(sample, type), by = c(sample_barcode = "sample")) %>%
  dplyr::left_join(tcga_surv, by = c(sample_barcode = "sample")) %>% 
  dplyr::mutate(sample_classification = factor(
    sample_classification, 
    c("No-fSCNA", "Linear", "Heavily-rearranged", "BFB", "Circular"))) %>% 
  dplyr::filter(!is.na(PFI.time), PFI.time <= 365 * 5)
fit = coxph(Surv(PFI.time, PFI) ~ sample_classification + strata(type), surv_df2)
forestmodel::forest_model(fit)

# PFI and Heavily-rearrangement may be missing points from NG paper


# ecGenes -----------------------------------------------------------------

load("modeling/data/train_data_extra_info_v3.RData")
colnames(data_gene_load)[2] = "ecGenes"

get_surv_data = function(rm_class = c("Linear", "Heavily-rearranged", "BFB")) {
  df = sample_tb %>% 
    dplyr::filter(startsWith(sample_barcode, "TCGA"), tumor_or_normal == "tumor") %>% 
    dplyr::select(sample_barcode, patient_barcode, lineage, sample_classification) %>% 
    dplyr::left_join(tcga_clinical %>%
                       dplyr::select(sample, type), by = c(sample_barcode = "sample")) %>%
    dplyr::left_join(tcga_surv, by = c(sample_barcode = "sample")) %>% 
    dplyr::mutate(sample_classification = factor(
      sample_classification, 
      c("No-fSCNA", "Linear", "Heavily-rearranged", "BFB", "Circular"))) %>% 
    dplyr::left_join(data_gene_load, data_gene_load, by = c("sample_barcode" = "sample")) %>% 
    dplyr::mutate(ecGenes = ifelse(is.na(ecGenes), 0, ecGenes/100), # normalize to unit per 100
                  ecStatus = ifelse(ecGenes > 0, 1L, 0L)) %>% 
    dplyr::filter(!sample_classification %in% rm_class) %>% dplyr::filter(!is.na(OS.time), OS.time <= 365 * 5)
  
  if (length(rm_class) == 1 && rm_class == "nonLinear") {
    df = df %>% 
      dplyr::mutate(nonlinear = ifelse(sample_classification %in%  c("No-fSCNA", "Linear"), 0L, 1L))
  }
  
  df
}


# Use subset of (No-fSCNA + Circular) and mostly No-fSCNA as reference
coxph(Surv(OS.time, OS) ~ ecGenes + strata(type), get_surv_data())
coxph(Surv(OS.time, OS) ~ ecStatus + strata(type), get_surv_data())


# Use subset of (No-fSCNA + Linear + Circular) and mostly No-fSCNA + Linear as reference
coxph(Surv(OS.time, OS) ~ ecGenes + strata(type), get_surv_data(c("Heavily-rearranged", "BFB"))) %>% 
  forestmodel::forest_model()
fit = coxph(Surv(OS.time, OS) ~ ecStatus + strata(type), get_surv_data(c("Heavily-rearranged", "BFB")))
forestmodel::forest_model(fit)

summary(fit)

df = get_surv_data(c("Heavily-rearranged", "BFB"))
table(class = df$sample_classification, ecStatus = df$ecStatus)

p = forestmodel::forest_model(fit, merge_models = TRUE)
p
ggplot2::ggsave("pancan-analysis/plots/tcga_ecStatus_cox_forest.pdf", p, width = 8, height = 3)

library(ggplot2)
p <- survminer::ggsurvplot(survfit(Surv(OS.time, OS) ~ ecStatus, 
                              data = get_surv_data(c("Heavily-rearranged", "BFB")) %>% 
                                dplyr::mutate(OS.time = OS.time / 365)),
                      pval = TRUE, fun = "pct", xlab = "Time (in years)",
                      palette = rev(c("red", "blue")),
                      legend = c(0.8, 0.9), 
                      legend.title = element_blank(), 
                      legend.labs = c("No-fSCNA/Linear", "ecDNA+"))
ggplot2::ggsave("pancan-analysis/plots/tcga_ecStatus_km_plot.pdf", p$plot, width = 6, height = 4)

# In survival aspcect, ecDNA exists or not is more important than the number of ecDNA genes
# 考虑肿瘤类型，相比于复杂重排和BFB，ecDNA组的生存不一定更差

