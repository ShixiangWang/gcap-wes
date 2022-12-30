# https://docs.google.com/spreadsheets/d/11VdviGiovgnSzTRb-0-Qv5SUfkXdoqSlAs4D9mFp1ao/edit#gid=0
# https://xenabrowser.net/datapages/?dataset=TCGASubtype.20170308.tsv&host=https%3A%2F%2Fpancanatlas.xenahubs.net&removeHub=https%3A%2F%2Fxena.treehouse.gi.ucsc.edu%3A443

library(UCSCXenaShiny)
library(gcaputils)

nrow(tcga_subtypes)
tcga_subtypes = data.table(tcga_subtypes)
tcga_subtypes = tcga_subtypes[, list(sample = sampleID, mol_subtype = Subtype_Selected, imu_subtype = Subtype_Immune_Model_Based)]

data_tcga = readRDS("data/tcga_snp_data.rds")
data_tcga$class = set_default_factor(data_tcga$class)
data_tcga = data_tcga[, list(sample, type, class, OS.time, OS, PFI.time, PFI)]
data = merge(data_tcga, tcga_subtypes, by = "sample")

table(is.na(data$OS.time), is.na(data$PFI.time))

# 按cancer type分析
library(survival)

## OS
type_rv = list()
for (t in unique(data$type)) {
  message("handling type ", t)
  type_rv[[t]] = list()
  type_mol = data[type == t & !is.na(mol_subtype) & !is.na(OS)]
  type_rv[[t]]$molecular_N = nrow(type_mol)
  type_rv[[t]]$molecular_table = table(type_mol$class)
  
  if (nrow(type_mol) > 20) {
    m1 = coxph(Surv(OS.time, OS) ~ mol_subtype, data = type_mol)
    m2 = coxph(Surv(OS.time, OS) ~ mol_subtype + class, data = type_mol)
    type_rv[[t]]$molecular_result = list(
      m1 = m1$concordance["concordance"],
      m2 = m2$concordance["concordance"],
      p = anova(m1, m2)$`P(>|Chi|)`[2]
    )
  }
  
  type_imu = data[type == t & !is.na(imu_subtype) & !is.na(OS)]
  type_rv[[t]]$immune_N = nrow(type_imu)
  type_rv[[t]]$immune_table = table(type_imu$class)
  
  if (nrow(type_imu) > 20) {
    m1 = coxph(Surv(OS.time, OS) ~ imu_subtype, data = type_imu)
    m2 = coxph(Surv(OS.time, OS) ~ imu_subtype + class, data = type_imu)
    type_rv[[t]]$immune_result = list(
      m1 = m1$concordance["concordance"],
      m2 = m2$concordance["concordance"],
      p = anova(m1, m2)$`P(>|Chi|)`[2]
    )
  }
}

type_rv_OS = type_rv

# PFS
type_rv = list()
for (t in unique(data$type)) {
  message("handling type ", t)
  type_rv[[t]] = list()
  type_mol = data[type == t & !is.na(mol_subtype) & !is.na(PFI)]
  type_rv[[t]]$molecular_N = nrow(type_mol)
  type_rv[[t]]$molecular_table = table(type_mol$class)
  
  if (nrow(type_mol) > 20) {
    m1 = coxph(Surv(PFI.time, PFI) ~ mol_subtype, data = type_mol)
    m2 = coxph(Surv(PFI.time, PFI) ~ mol_subtype + class, data = type_mol)
    type_rv[[t]]$molecular_result = list(
      m1 = m1$concordance["concordance"],
      m2 = m2$concordance["concordance"],
      p = anova(m1, m2)$`P(>|Chi|)`[2]
    )
  }
  
  type_imu = data[type == t & !is.na(imu_subtype) & !is.na(PFI)]
  type_rv[[t]]$immune_N = nrow(type_imu)
  type_rv[[t]]$immune_table = table(type_imu$class)
  
  if (nrow(type_imu) > 20) {
    m1 = coxph(Surv(PFI.time, PFI) ~ imu_subtype, data = type_imu)
    m2 = coxph(Surv(PFI.time, PFI) ~ imu_subtype + class, data = type_imu)
    type_rv[[t]]$immune_result = list(
      m1 = m1$concordance["concordance"],
      m2 = m2$concordance["concordance"],
      p = anova(m1, m2)$`P(>|Chi|)`[2]
    )
  }
}

type_rv_PFI = type_rv

sapply(type_rv_PFI, function(x) x$immune_result$p)

save(type_rv_OS, type_rv_PFI, file = "data/TCGA_subtypes_analysis.RData")

rv_os_mol = list(
  m2 = unlist(sapply(type_rv_OS, function(x) x$molecular_result$m2 - x$molecular_result$m1)),
  p = unlist(sapply(type_rv_OS, function(x) x$molecular_result$p))
)
rv_os_imu = list(
  m2 = unlist(sapply(type_rv_OS, function(x) x$immune_result$m2 - x$immune_result$m1)),
  p = unlist(sapply(type_rv_OS, function(x) x$immune_result$p))
)

rv_pfs_mol = list(
  m2 = unlist(sapply(type_rv_PFI, function(x) x$molecular_result$m2 - x$molecular_result$m1)),
  p = unlist(sapply(type_rv_PFI, function(x) x$molecular_result$p))
)
rv_pfs_imu = list(
  m2 = unlist(sapply(type_rv_PFI, function(x) x$immune_result$m2 - x$immune_result$m1)),
  p = unlist(sapply(type_rv_PFI, function(x) x$immune_result$p))
)

type_result = rbind(
  data.table(g1 = "OS", g2 = "Molecular", type = names(rv_os_mol$p), c.index.increase = as.numeric(rv_os_mol$m2), p = as.numeric(rv_os_mol$p)),
  data.table(g1 = "PFS", g2 = "Molecular", type = names(rv_pfs_mol$p), c.index.increase = as.numeric(rv_pfs_mol$m2), p = as.numeric(rv_pfs_mol$p)),
  data.table(g1 = "OS", g2 = "Immune", type = names(rv_os_imu$p), c.index.increase = as.numeric(rv_os_imu$m2), p = as.numeric(rv_os_imu$p)),
  data.table(g1 = "PFS", g2 = "Immune", type = names(rv_pfs_imu$p), c.index.increase = as.numeric(rv_pfs_imu$m2), p = as.numeric(rv_pfs_imu$p))
)

saveRDS(type_result, file = "data/TCGA_subtype_analysis_tidy_result.rds")

library(ggplot2)
ggplot(type_result %>% 
         dplyr::mutate(group = paste(g2, g1, sep = " subtype - "),
                       label = dplyr::case_when(
                         p < 0.001 ~ "***",
                         p < 0.01 ~ "**",
                         p < 0.05 ~ "*",
                         TRUE ~ ""
                       )), 
       aes(x = type, y = group)) +
  geom_tile(aes(fill = c.index.increase)) +
  scale_fill_binned(type = "viridis") +
  geom_text(aes(label = label), color = "white") +
  labs(x = NULL, y = NULL) +
  theme_base() +
  #cowplot::theme_cowplot() +
  theme(legend.position = "top", legend.text = element_text(size = 8)) +
  ggpubr::rotate_x_text() -> p
p
ggsave("plots/TCGA_subtype_add_fCNA_analysis_heatmap.pdf", p, width = 9, height = 4)

# 做一个带CIN的肿瘤的分析，即消化道肿瘤
data2 = data[type %in% c("ESCA", "STAD", "COAD", "READ")]
table(data2$mol_subtype)
# GI.ESCC?
# same class from tcgabiolinks
# subtypes <- PanCancerAtlas_subtypes()
data2[, mol_subtype := ifelse(mol_subtype %in% c("GI.HM-indel", "GI.HM-SNV"), "GI.HM", mol_subtype)]
table(data2$mol_subtype)

devtools::load_all("~/proj/gcaputils/")
#debug(gcap.plotKMcurve)
gcap.plotKMcurve(data2[, list(sample, mol_subtype, OS.time, OS)], palette = RColorBrewer::brewer.pal(5, "Set1"))
gcap.plotKMcurve(data2[, list(sample, mol_subtype, PFI.time, PFI)], palette = RColorBrewer::brewer.pal(5, "Set1"))


m1 = coxph(Surv(OS.time, OS) ~ mol_subtype, data = data2)
m2 = coxph(Surv(OS.time, OS) ~ mol_subtype + class, data = data2)
list(
  m1 = m1$concordance["concordance"],
  m2 = m2$concordance["concordance"],
  p = anova(m1, m2)$`P(>|Chi|)`[2]
) -> gi_1

m1 = coxph(Surv(PFI.time, PFI) ~ mol_subtype, data = data2)
m2 = coxph(Surv(PFI.time, PFI) ~ mol_subtype + class, data = data2)
list(
  m1 = m1$concordance["concordance"],
  m2 = m2$concordance["concordance"],
  p = anova(m1, m2)$`P(>|Chi|)`[2]
) -> gi_2

m1 = coxph(Surv(OS.time, OS) ~ imu_subtype, data = data2)
m2 = coxph(Surv(OS.time, OS) ~ imu_subtype + class, data = data2)
list(
  m1 = m1$concordance["concordance"],
  m2 = m2$concordance["concordance"],
  p = anova(m1, m2)$`P(>|Chi|)`[2]
) -> gi_3

m1 = coxph(Surv(PFI.time, PFI) ~ imu_subtype, data = data2)
m2 = coxph(Surv(PFI.time, PFI) ~ imu_subtype + class, data = data2)
list(
  m1 = m1$concordance["concordance"],
  m2 = m2$concordance["concordance"],
  p = anova(m1, m2)$`P(>|Chi|)`[2]
) -> gi_4

gi_rv = rbind(
  cbind(g1 = "molecular", g2 = "OS", as.data.table(gi_1)),
  cbind(g1 = "molecular", g2 = "PFS", as.data.table(gi_2)),
  cbind(g1 = "immune", g2 = "OS", as.data.table(gi_3)),
  cbind(g1 = "immune", g2 = "PFS", as.data.table(gi_4))
)

gi_rv
saveRDS(gi_rv, file = "data/TCGA_GI_subtype_analysis_result.rds")

gi2 = melt(gi_rv, id.vars = c("g1", "g2"), measure.vars = c("m1", "m2"))

# ggplot(gi2 %>% dplyr::mutate(g1 = paste(g1, "subtype")), 
#        aes(x = variable, y = value, fill = variable)) +
#   geom_bar(stat = "identity", width = 0.5) +
#   facet_grid(g1 ~ g2) +
#   cowplot::theme_cowplot() +
#   ylim(0, 0.7) + labs(x = NULL, y = "C index") +
#   theme(legend.position = "none") -> p
# p
ggplot(gi2[g1 == "molecular"], 
       aes(x = variable, y = value, fill = variable)) +
  geom_bar(stat = "identity", width = 0.5) +
  facet_grid(. ~ g2) +
  cowplot::theme_cowplot() +
  ylim(0, 0.7) + labs(x = NULL, y = "C index") +
  theme(legend.position = "none") -> p
p
ggsave("plots/TCGA_GI_c_index.pdf", p, width = 4, height = 3)

data2[, combined := paste0(mol_subtype, "/", class)]
table(data2$combined)
data2[, subtype := fcase(
  mol_subtype == "GI.HM", "HM",
  mol_subtype == "GI.EBV", "EBV",
  combined %in% c("GI.CIN/nofocal"), "CIN-Mild",
  combined %in% c("GI.GS/nofocal", "GI.ESCC/nofocal"), "CN-Quiet",
  class == "noncircular", "Non-Circ",
  class == "circular", "Circ",
  default = "CN-Quiet"
)]
new_lvls = c("HM", "CIN-Mild", "CN-Quiet", "Non-Circ", "EBV", "Circ")
data2$subtype = factor(data2$subtype, levels = new_lvls)
new_pal_cols = c("grey", "black", "cyan", "blue", "magenta", "red")
# RColorBrewer::brewer.pal(6, "Set1")

gcap.plotKMcurve(data2[, list(sample, subtype, OS.time, OS)], palette = new_pal_cols)
gcap.plotKMcurve(data2[, list(sample, subtype, PFI.time, PFI)], palette = new_pal_cols)

# 其实这个分析只需要关注 GS/CIN 再分型fCNA的效果就可以
p1 = gcap.plotKMcurve(data2[mol_subtype == "GI.GS", list(sample, class, OS.time, OS)])
p2 = gcap.plotKMcurve(data2[mol_subtype == "GI.GS", list(sample, class, PFI.time, PFI)])
p3 = gcap.plotKMcurve(data2[mol_subtype == "GI.CIN", list(sample, class, OS.time, OS)])
p4 = gcap.plotKMcurve(data2[mol_subtype == "GI.CIN", list(sample, class, PFI.time, PFI)])
p5 = gcap.plotKMcurve(data2[mol_subtype %in% c("GI.CIN", "GI.GS"), list(sample, class, OS.time, OS)])
p6 = gcap.plotKMcurve(data2[mol_subtype %in% c("GI.CIN", "GI.GS"), list(sample, class, PFI.time, PFI)])


pdf("plots/TCGA_GI_GS_fCNA_OS.pdf", width = 6, height = 6, onefile = FALSE)
p1
dev.off()

pdf("plots/TCGA_GI_GS_fCNA_PFS.pdf", width = 6, height = 6, onefile = FALSE)
p2
dev.off()

pdf("plots/TCGA_GI_CIN_fCNA_OS.pdf", width = 6, height = 6, onefile = FALSE)
p3
dev.off()

pdf("plots/TCGA_GI_CIN_fCNA_PFS.pdf", width = 6, height = 6, onefile = FALSE)
p4
dev.off()

pdf("plots/TCGA_GI_GS_plus_CIN_fCNA_OS.pdf", width = 6, height = 6, onefile = FALSE)
p5
dev.off()

pdf("plots/TCGA_GI_GS_plus_CIN_fCNA_PFS.pdf", width = 6, height = 6, onefile = FALSE)
p6
dev.off()


