library(gcaputils)

data_tcga = readRDS("data/tcga_snp_data.rds")
tcga_crc = data_tcga[type %in% c("COAD", "READ")]
tcga_crc$class = set_default_factor(tcga_crc$class)

tcga_crc$class %>% table()

p1 = gcap.plotKMcurve(tcga_crc[, .(sample, class, OS.time, OS)])
p2 = gcap.plotKMcurve(tcga_crc[, .(sample, class, PFI.time, PFI)])
p3 = gcap.plotKMcurve(tcga_crc[, .(sample, class, DFI.time, DFI)])


pdf("plots/TCGA_CRC_for_fCNA_class_OS.pdf", width = 7, height = 7, onefile = FALSE)
print(p1)
dev.off()

pdf("plots/TCGA_CRC_for_fCNA_class_PFS.pdf", width = 7, height = 7, onefile = FALSE)
print(p2)
dev.off()

pdf("plots/TCGA_CRC_for_fCNA_class_DFS.pdf", width = 7, height = 7, onefile = FALSE)
print(p3)
dev.off()
