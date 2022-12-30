# Using SNP6 microarray data, copy number profiles were generated for 9873 cancers and
# matching germline DNA of 33 different types from TCGA using ASCAT with a segmentation
# penalty of 70. (Steele, Nature, 2022)
#
# 9699 cancers successfully called

setwd(file.path(PROJ_DIR, "preprocessing-and-EDA/"))

ascat = readxl::read_excel("/data3/wsx_data/Nature2022_CNS/Nature2022_ASCAT.xlsx", sheet = 2) %>% data.table()
info = readxl::read_excel("/data3/wsx_data/Nature2022_CNS/Nature2022_ASCAT.xlsx", sheet = 3) %>% data.table()

ascat2 = merge(ascat, info, by.x = "sample", by.y = "name")
ascat2[, barcode := substr(barcodeTumour, 1, 15)]
ascat2[, total_cn := nMajor + nMinor]
ascat2 = ascat2[, -c("sample", "patient", "CNclass", "barcodeTumour", "barcodeNormal", "solution", "frac_homo",
            "tumour_mapd", "normal_mapd", "segDiff", "pass", "rep", "nMajor")]


data.table::setnames(ascat2, c("startpos", "endpos", "nMinor", "barcode"),
                     c("start", "end", "minor_cn", "sample"))

data.table::setcolorder(ascat2, c("sample", "chr", "start", "end", "total_cn", "minor_cn"))

uniqLen(ascat2$sample)
# 9699

saveRDS(ascat2, file = "data/tcga_ascat.rds")

